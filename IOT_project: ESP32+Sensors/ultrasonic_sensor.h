// ------------- UltarSonic sensors START -------------------
// Include NewPing Library
#define trigPin 21
#define echoPin 39



//#define MAX_DISTANCE 400.0 //cm
//#define MIN_DISTANCE 2.0 //cm
// use it to check the level of the garbage
float LOCAL_STATE = 0;

// save the last state we uploded to the cloud
float SERVER_STATE = 0; 



volatile float bin_depth_cm = -1;
float bin_level_cm = -1;
float max_distance = 400;
float THRASHOLD = 400;

// update new max distance
void update_max_distance_and_thrashhold(){
  max_distance = min( float(400), float(bin_depth_cm));
  THRASHOLD = max_distance*0.1;
}

void update_last_state(){
  SERVER_STATE = LOCAL_STATE;
}

// func to check if garbage level changed at a quntity justifying updating the server
bool check_if_update(){

  if (bin_depth_cm<0) return false;
  float sum = fabs(SERVER_STATE - LOCAL_STATE);
  // right case at '||' is where we are at the last THRASHOLD level, 90%< so update each change of 5 cm
  if (sum > THRASHOLD || ((LOCAL_STATE <=THRASHOLD) && sum > 5)){
    return true; 
  }
  return false;
}


bool read_ultasonic_sensor(){
  digitalWrite(trigPin, LOW); // set to 0 
  delayMicroseconds(2); // wait 2 us

  digitalWrite(trigPin, HIGH); // set to 1
  //delay(10); // wait 10 us
  delayMicroseconds(10); 
  digitalWrite(trigPin, LOW); // set to low

// Reads a pulse (either HIGH or LOW) on a pin.
// if value is HIGH, pulseIn() waits for the pin to go from LOW to HIGH, starts timing, then waits for the pin to go LOW and stops timing.
// Returns the length of the pulse in microseconds or gives up and returns 0 if no complete pulse was received within the timeout.
  float duration = pulseIn(echoPin, HIGH);
  
  // Determine distance from duration
  // Use 346 metres per second as speed of sound at 25 c
  // divide 346 by 10^4 since we want cm / microsec
  // divide by 2 for the distance of sensor from object
  float distance = (duration / 2) * 0.0346;

  //make sure distance is smaller\ equal to bin_depth && 400 cm, and greater then 2 cm
  if (distance > max_distance  || distance <= 2.0) {
    Serial.println("Out of range");
    return false;
  }

  else {
    Serial.print("> ultrasonic sensor measurement is: ");
    Serial.print(distance); // print the distance
    Serial.println(" cm"); 
  }

  LOCAL_STATE = distance;
  bin_level_cm = bin_depth_cm - LOCAL_STATE;
  return true;
}

// ---------------- will delete -------------------


bool read_ultasonic_sensor_iter(){
  bool check = false;
  for (int i=0; i<5; ++i){
    check = read_ultasonic_sensor();
    if (check) {
      break;
    }

  }
  return check;
}

