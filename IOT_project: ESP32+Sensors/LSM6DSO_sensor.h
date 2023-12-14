// ------------- LSM6DSO START ---------------

#include "SparkFunLSM6DSO.h"
#include "Wire.h"
//#include "SPI.h"

LSM6DSO myIMU; //Default constructor is I2C, addr 0x6B
uint8_t reg;
void initalize_LSM(){
  delay(500);  
  Wire.begin();
  delay(10);
  if( myIMU.begin() )
    Serial.println("Ready.");
  else { 
    Serial.println("Could not connect to IMU.");
    Serial.println("Freezing");
  }

  if( myIMU.initialize(BASIC_SETTINGS) )
    Serial.println("Loaded Settings.");
};


float calc_temperature(){
  //Serial.print("\nThermometer:\n");
  //Serial.print(" Degrees C = ");
  // read Celius tempreture
  return myIMU.readTempC();
  //Serial.println(temperature, 3);
};

void Wake_Up_Mode(){

  myIMU.writeRegister(0x5C, 0x10);  // No duration time and WAKE_THS_W = 1, 1 LSB = FS_XL / 2^8
  delay(10);
  myIMU.writeRegister(0x5B, 0x08);  // Set wake-up threshold. 
  delay(10);
  myIMU.writeRegister(0x56, 0x50);  // Enable interrupts and apply slope filter; latch mode disabled <=> LIR bit=0
  delay(10);
  myIMU.writeRegister(0x58,  0x80);   // Enable interrupt function 
  delay(10);
  myIMU.writeRegister(0x10,  0x70);   // Turn on the accelerometer, ODR_XL = 833 Hz, FS_XL = ±2 g
  
  delay(50);                          //this delay time after turning on the accelerometer and before enable the signal to INT1 pin solve the edge case of first sample 

  myIMU.writeRegister( 0x10, 0x10 );  // change to slower rate, ODR_XL = 12.5 Hz
  delay(10);
  myIMU.writeRegister(0x5E, 0x20);  // Wake-up interrupt driven to INT1 pin
}


void initalize_LSM_with_interrupt() {

  
  Serial.begin(115200);
  delay(500); 
  Wire.begin();
  delay(10);
  if( myIMU.begin() )
    Serial.println("Ready.");
  else { 
    Serial.println("Could not connect to IMU.");
    Serial.println("Freezing");
  }
  delay(10);
  Wake_Up_Mode();
  delay(10);
}




// ===================== not gonna use it =======================================================================================
void relitive_tilt(){
 
  myIMU.writeRegister( 0x01, 0x80 ); // Enable access to embedded functions registers // Enable significant motion detection
  delay(10);
  myIMU.writeRegister( 0x04, 0x10 ); // Enable significant motion detection 
  delay(10);
  myIMU.writeRegister( 0x0A, 0x10 ); // Significant motion interrupt driven to INT1 pin 
  delay(10);
  myIMU.writeRegister( 0x17, 0x80 ); // Enable latched mode for embedded functions  
  delay(10);
  myIMU.writeRegister( 0x01, 0x00 ); // Disable access to embedded functions registers 
  delay(10);
  myIMU.writeRegister( 0x5E, 0x02 ); // Enable embedded functions interrupt routing
  delay(10);
  myIMU.writeRegister( 0x10, 0x20 ); //  CTRL1_XL = 0x20.Turn on the accelerometer => ODR_XL = 26 Hz and FS_XL = ±2 g
  delay(10);

}
void reset_latch(){
  digitalRead(GPIO_NUM_15);
  return;
}
// ===============================================================================================================================
