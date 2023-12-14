// Copyright (c) Microsoft Corporation. All rights reserved.
// SPDX-License-Identifier: MIT

/*
 * This is an Arduino-based Azure IoT Hub sample for ESPRESSIF ESP32 boards.
 * It uses our Azure Embedded SDK for C to help interact with Azure IoT.
 * For reference, please visit https://github.com/azure/azure-sdk-for-c.
 *
 * To connect and work with Azure IoT Hub you need an MQTT client, connecting, subscribing
 * and publishing to specific topics to use the messaging features of the hub.
 * Our azure-sdk-for-c is an MQTT client support library, helping composing and parsing the
 * MQTT topic names and messages exchanged with the Azure IoT Hub.
 *
 * This sample performs the following tasks:
 * - Synchronize the device clock with a NTP server;
 * - Initialize our "az_iot_hub_client" (struct for data, part of our azure-sdk-for-c);
 * - Initialize the MQTT client (here we use ESPRESSIF's esp_mqtt_client, which also handle the tcp
 * connection and TLS);
 * - Connect the MQTT client (using server-certificate validation, SAS-tokens for client
 * authentication);
 * - Periodically send telemetry data to the Azure IoT Hub.
 *
 * To properly connect to your Azure IoT Hub, please fill the information in the `iot_configs.h`
 * file.
 */
#include <ArduinoJson.h>
#include "ultrasonic_sensor.h"
#include "LSM6DSO_sensor.h"
// C99 libraries
#include <cstdlib>
#include <string.h>
#include <time.h>

// Libraries for MQTT client and WiFi connection
#include <WiFi.h>
#include <mqtt_client.h>

// Azure IoT SDK for C includes
#include <az_core.h>
#include <az_iot.h>
#include <azure_ca.h>

// Additional sample headers
#include "AzIoTSasToken.h"
#include "SerialLogger.h"
#include "iot_configs.h"

// When developing for your own Arduino-based platform,
// please follow the format '(ard;<platform>)'.
#define AZURE_SDK_CLIENT_USER_AGENT "c%2F" AZ_SDK_VERSION_STRING "(ard;esp32)"

// Utility macros and defines
#define sizeofarray(a) (sizeof(a) / sizeof(a[0]))
#define NTP_SERVERS "pool.ntp.org", "time.nist.gov"
#define MQTT_QOS1 1
#define DO_NOT_RETAIN_MSG 0
#define SAS_TOKEN_DURATION_IN_MINUTES 60
#define UNIX_TIME_NOV_13_2017 1510592825

#define PST_TIME_ZONE -8
#define PST_TIME_ZONE_DAYLIGHT_SAVINGS_DIFF 1

#define GMT_OFFSET_SECS (PST_TIME_ZONE * 3600)
#define GMT_OFFSET_SECS_DST ((PST_TIME_ZONE + PST_TIME_ZONE_DAYLIGHT_SAVINGS_DIFF) * 3600)

// Translate iot_configs.h defines into variables used by the sample
static const char* ssid = IOT_CONFIG_WIFI_SSID;
static const char* password = IOT_CONFIG_WIFI_PASSWORD;
static const char* host = IOT_CONFIG_IOTHUB_FQDN;
static const char* mqtt_broker_uri = "mqtts://" IOT_CONFIG_IOTHUB_FQDN;
static const char* device_id = IOT_CONFIG_DEVICE_ID;
static const int mqtt_port = AZ_IOT_DEFAULT_MQTT_CONNECT_PORT;

// Memory allocated for the sample's variables and structures.
static esp_mqtt_client_handle_t mqtt_client;
static az_iot_hub_client client;

static char mqtt_client_id[128];
static char mqtt_username[128];
static char mqtt_password[200];
static uint8_t sas_signature_buffer[256];
static unsigned long next_telemetry_send_time_ms = 0;
//static char telemetry_topic[128];
//static uint8_t telemetry_payload[100];
static uint32_t telemetry_send_count = 0;

#define INCOMING_DATA_BUFFER_SIZE 128
static char incoming_data[INCOMING_DATA_BUFFER_SIZE];

// Auxiliary functions
#ifndef IOT_CONFIG_USE_X509_CERT
static AzIoTSasToken sasToken(
    &client,
    AZ_SPAN_FROM_STR(IOT_CONFIG_DEVICE_KEY),
    AZ_SPAN_FROM_BUFFER(sas_signature_buffer),
    AZ_SPAN_FROM_BUFFER(mqtt_password));
#endif // IOT_CONFIG_USE_X509_CERT

static void connectToWiFi()
{
  //Logger.Info("Connecting to WIFI SSID " + String(ssid));
    Serial.println("> Connecting to WIFI...");

  WiFi.mode(WIFI_STA);
  WiFi.disconnect();
  delay(100);
  WiFi.begin(ssid, password);
  while (WiFi.status() != WL_CONNECTED)
  {
    delay(500);
    Serial.print(".");
  }

  Serial.println("");
  Serial.println("> WiFi connected");
  //Logger.Info("WiFi connected, IP address: " + WiFi.localIP().toString());
}


// 
static int connectToWiFi_2()
{
  //Logger.Info("Connecting to WIFI SSID " + String(ssid));
    Serial.println("> Connecting to WIFI...");

  WiFi.mode(WIFI_STA);
  WiFi.disconnect();
  delay(100);
  WiFi.begin(ssid, password);
  for(int k=0; k<180; k++){
    delay(500);
    if(WiFi.status() == WL_CONNECTED){
      Serial.println("");
      Serial.println("> WiFi connected");
      return 0;
    }
    
  }
      Serial.println("");
      Serial.println("> unable to conncet to WiFi");
      return 1;
}


static void initializeTime()
{
  //Logger.Info("Setting time using SNTP");

  configTime(GMT_OFFSET_SECS, GMT_OFFSET_SECS_DST, NTP_SERVERS);
  time_t now = time(NULL);
  while (now < UNIX_TIME_NOV_13_2017)
  {
    delay(500);
    Serial.print(".");
    now = time(nullptr);
  }
  Serial.println("");
  //Logger.Info("Time initialized!");
  Serial.println("> Time initialized");
}

void receivedCallback(char* topic, byte* payload, unsigned int length)
{
  Logger.Info("Received [");
  Logger.Info(topic);
  Logger.Info("]: ");
  for (int i = 0; i < length; i++)
  {
    Serial.print((char)payload[i]);
  }
  Serial.println("");
}

static esp_err_t mqtt_event_handler(esp_mqtt_event_handle_t event)
{
  switch (event->event_id)
  {
    int i, r;

    case MQTT_EVENT_ERROR:
      Logger.Info("MQTT event MQTT_EVENT_ERROR");
      break;
    case MQTT_EVENT_CONNECTED:
      //Logger.Info("MQTT event MQTT_EVENT_CONNECTED");

      r = esp_mqtt_client_subscribe(mqtt_client, AZ_IOT_HUB_CLIENT_C2D_SUBSCRIBE_TOPIC, 1);
      if (r == -1)
      {
        Logger.Error("Could not subscribe for cloud-to-device messages.");
      }
      else
      {
        //Logger.Info("Subscribed for cloud-to-device messages; message id:" + String(r));
      }

      break;
    case MQTT_EVENT_DISCONNECTED:
      //Logger.Info("MQTT event MQTT_EVENT_DISCONNECTED");
      break;
    case MQTT_EVENT_SUBSCRIBED:
      //Logger.Info("MQTT event MQTT_EVENT_SUBSCRIBED");
      break;
    case MQTT_EVENT_UNSUBSCRIBED:
      //Logger.Info("MQTT event MQTT_EVENT_UNSUBSCRIBED");
      break;
    case MQTT_EVENT_PUBLISHED:
      //Logger.Info("MQTT event MQTT_EVENT_PUBLISHED");
      break;
    case MQTT_EVENT_DATA:
      //Logger.Info("MQTT event MQTT_EVENT_DATA");


      //for (i = 0; i < (INCOMING_DATA_BUFFER_SIZE - 1) && i < event->topic_len; i++)
      //{
      //  incoming_data[i] = event->topic[i];
      //}
      //incoming_data[i] = '\0';
      //Logger.Info("Topic: " + String(incoming_data));
     

      for (i = 0; i < (INCOMING_DATA_BUFFER_SIZE - 1) && i < event->data_len; i++)
      {
        incoming_data[i] = event->data[i];
      }
      incoming_data[i] = '\0';
      //String payloadMessage = String(incoming_data);
      //Logger.Info("Data: " + String(incoming_data));
      message_handler(incoming_data);
      
      break;
    case MQTT_EVENT_BEFORE_CONNECT:
      //Logger.Info("MQTT event MQTT_EVENT_BEFORE_CONNECT");
      break;
    default:
      Logger.Error("MQTT event UNKNOWN");
      break;
  }   
  return ESP_OK;
}

static void initializeIoTHubClient()
{
  az_iot_hub_client_options options = az_iot_hub_client_options_default();
  options.user_agent = AZ_SPAN_FROM_STR(AZURE_SDK_CLIENT_USER_AGENT);

  if (az_result_failed(az_iot_hub_client_init(
          &client,
          az_span_create((uint8_t*)host, strlen(host)),
          az_span_create((uint8_t*)device_id, strlen(device_id)),
          &options)))
  {
    Logger.Error("Failed initializing Azure IoT Hub client");
    return;
  }

  size_t client_id_length;
  if (az_result_failed(az_iot_hub_client_get_client_id(
          &client, mqtt_client_id, sizeof(mqtt_client_id) - 1, &client_id_length)))
  {
    Logger.Error("Failed getting client id");
    return;
  }

  if (az_result_failed(az_iot_hub_client_get_user_name(
          &client, mqtt_username, sizeofarray(mqtt_username), NULL)))
  {
    Logger.Error("Failed to get MQTT clientId, return code");
    return;
  }

  //Logger.Info("Client ID: " + String(mqtt_client_id));
  //Logger.Info("Username: " + String(mqtt_username));
}

static int initializeMqttClient()
{
#ifndef IOT_CONFIG_USE_X509_CERT
  if (sasToken.Generate(SAS_TOKEN_DURATION_IN_MINUTES) != 0)
  {
    Logger.Error("Failed generating SAS token");
    return 1;
  }
#endif

  esp_mqtt_client_config_t mqtt_config;
  memset(&mqtt_config, 0, sizeof(mqtt_config));
  mqtt_config.uri = mqtt_broker_uri;
  mqtt_config.port = mqtt_port;
  mqtt_config.client_id = mqtt_client_id;
  mqtt_config.username = mqtt_username;

#ifdef IOT_CONFIG_USE_X509_CERT
  Logger.Info("MQTT client using X509 Certificate authentication");
  mqtt_config.client_cert_pem = IOT_CONFIG_DEVICE_CERT;
  mqtt_config.client_key_pem = IOT_CONFIG_DEVICE_CERT_PRIVATE_KEY;
#else // Using SAS key
  mqtt_config.password = (const char*)az_span_ptr(sasToken.Get());
#endif

  mqtt_config.keepalive = 30;
  mqtt_config.disable_clean_session = 0;
  mqtt_config.disable_auto_reconnect = false;
  mqtt_config.event_handle = mqtt_event_handler;
  mqtt_config.user_context = NULL;
  mqtt_config.cert_pem = (const char*)ca_pem;

  mqtt_client = esp_mqtt_client_init(&mqtt_config);

  if (mqtt_client == NULL)
  {
    Logger.Error("Failed creating mqtt client");
    return 1;
  }

  esp_err_t start_result = esp_mqtt_client_start(mqtt_client);

  if (start_result != ESP_OK)
  {
    Logger.Error("Could not start mqtt client; error code:" + start_result);
    return 1;
  }
  else
  {
    //Logger.Info("MQTT client started");
    return 0;
  }
}

/*
 * @brief           Gets the number of seconds since UNIX epoch until now.
 * @return uint32_t Number of seconds.
 */
static uint32_t getEpochTimeInSecs() { return (uint32_t)time(NULL); }

static void establishConnection()
{
  connectToWiFi();
  initializeTime();
  initializeIoTHubClient();
  (void)initializeMqttClient();
}


// ------------------------ add ---------------------------

void message_handler(char* incoming_data){
    String payload = String(incoming_data);
    Serial.println("> Message received: " +  payload); // my addition
    // convert message to json
    StaticJsonDocument<200> doc;
    DeserializationError error = deserializeJson(doc, payload);
      
    // print the json document
    //serializeJsonPretty(doc, Serial);
    if (doc["c2d_message_type"] == "update_bin_depth"){
       // output to the console to wait for the counters to be updated
      Serial.println("> Waiting for the bin depth to be updated... Please wait."); 
      // check that the message contains the counter id and counter value
      if  (!doc.containsKey("bin_depth_cm"))
      {
        Serial.println("> Message does not contain the bin depth...");
        //print_menu();
        return;
      }
        // Check if the "age" field is a float
      if (!doc["bin_depth_cm"].is<float>()) {
          Serial.println("> Bin depth is not a float...");
          //print_menu();
          return;
      }
      bin_depth_cm = doc["bin_depth_cm"];
      bin_level_cm = bin_depth_cm;
      // output to the console that the counters have been updated
      Serial.println("> Bin depth has been updated.");
    }
    else
    {
    Serial.println("> Message not recognized...");
    }    
    //print_menu();
    return;
}


int get_bin_height(){
    if (WiFi.status() != WL_CONNECTED)
    {
    connectToWiFi();
    }
    #ifndef IOT_CONFIG_USE_X509_CERT
    else if (sasToken.IsExpired())
    {
    Logger.Info("SAS token expired; reconnecting with a new one.");
    (void)esp_mqtt_client_destroy(mqtt_client);
    initializeMqttClient();
    }
    #endif
    char telemetry_topic[256];
    // construct device 2 cloud message
     StaticJsonDocument<200> doc;
    doc["bin_id"] = IOT_CONFIG_DEVICE_ID;
    doc["d2c_message_type"] = "get_bin_depth";


    char payload[256];
    int serializedSize = serializeJson(doc, payload);
      // publish message
    
    if (az_result_failed(az_iot_hub_client_telemetry_get_publish_topic(
            &client, NULL, telemetry_topic, sizeof(telemetry_topic), NULL)))
    {
        Serial.println("> Failed az_iot_hub_client_telemetry_get_publish_topic");
        return 2;
    }

    if (esp_mqtt_client_publish(
          mqtt_client,
          telemetry_topic,
           (const char*)payload,
          serializedSize,
          MQTT_QOS1,
          DO_NOT_RETAIN_MSG)
      == 0)
  {
    Logger.Error("Failed publishing");
    return 3;
  }
  else
  {
    Logger.Info("Message published successfully");
  }
    return 0;
}


int update_bin_level(){
    if (WiFi.status() != WL_CONNECTED)
    {
    int ret = connectToWiFi_2();
    if(ret){
      return 1;
    }
    }
  #ifndef IOT_CONFIG_USE_X509_CERT
  else if (sasToken.IsExpired())
  {
    Logger.Info("SAS token expired; reconnecting with a new one.");
    (void)esp_mqtt_client_destroy(mqtt_client);
    initializeMqttClient();
  }
  #endif
  
    char telemetry_topic[256];
    // construct device 2 cloud message
     StaticJsonDocument<200> doc;
    doc["bin_id"] = device_id;
    doc["bin_level_cm"] = bin_level_cm;
    doc["d2c_message_type"] = "update_bin_level";
    char payload[256];
    int serializedSize = serializeJson(doc, payload);
      // publish message
    
    if (az_result_failed(az_iot_hub_client_telemetry_get_publish_topic(
            &client, NULL, telemetry_topic, sizeof(telemetry_topic), NULL)))
    {
        Serial.println("> Failed az_iot_hub_client_telemetry_get_publish_topic");
        return 2;
    }

    if (esp_mqtt_client_publish(
          mqtt_client,
          telemetry_topic,
           (const char*)payload,
          serializedSize,
          MQTT_QOS1,
          DO_NOT_RETAIN_MSG)
      == 0)
  {
    Logger.Error("Failed publishing");
    return 3;
  }
  else
  {
    Logger.Info("Message published successfully");
  }
    return 0;
}


void setup() {
  delay(250);
  Serial.begin(115200);
  delay(250);
  Serial.println("Initializing device, please wait...");
  delay(50);
  // set input/output pins
  pinMode(GPIO_NUM_15, INPUT);
  pinMode(trigPin, OUTPUT); 
  pinMode(echoPin, INPUT);
  delay(50);

  // set pin 15 as interupt pin and check validity
  esp_err_t t = esp_sleep_enable_ext0_wakeup(GPIO_NUM_15, 1); 
  if(t!= ESP_OK){
    Serial.print("bad pin number: ");
    Serial.println(t);
  }
  delay(50);

  // setup sensor
  initalize_LSM_with_interrupt();
  delay(100);
  establishConnection();
  delay(10);
  get_bin_height();
  //read_ultasonic_sensor(); // calc distance
  //update_bin_level();
  
  //delay 10 sec to recieve messages
  delay(10000); 
  update_max_distance_and_thrashhold();
}


void loop()
{
  delay(10);
  Serial.println("going to sleep, wake me up if someone open the can...");
  delay(1000);  //delay 1 sec
  esp_light_sleep_start();
  Serial.println("I'm awake! checking can level...");
  delay(7000); //delay 7 sec. at real life we might change it to a greater value

  if( read_ultasonic_sensor() ){ // if result is out of range return false and go to sleep, else check if update is needed
      if ( check_if_update() ){ // check if there is a change at bin level greater then threshold
        if( update_bin_level() ){ // try to update the server with the new data. 0 at seccues else int greater then 0.
         Serial.println("unable to update server with the bin level");
      } 
    }
  } 
}
