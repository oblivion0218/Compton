#include <Servo.h>
#define SERVO1 23 // Alto
#define SERVO2 29 // Basso

Servo myservo1;  // create servo object to control a servo
Servo myservo2;  // create servo object to control a servo

void setup() {
// Servo motore con potenziometro
  myservo1.attach(SERVO1);  // attaches the servo on pin 9 to the servo object
  myservo2.attach(SERVO2);  // attaches the servo on pin 9 to the servo object
}

void loop() {
// // Servo motore con potenziometro
  // int input = analogRead(POTENZIOMETRO);            // reads the value of the potentiometer (value between 0 and 1023)
  // int angolo = map(input, 0, 1023, 0, 180);     // scale it to use it with the servo (value between 0 and 180)

  myservo1.write(23 + 0);     //23 x 90        // sets the servo position according to the scaled value
  myservo2.write(59 + 90);      //59 x 90            // sets the servo position according to the scaled value

  // for(int i = 0; i <= 90; i += 10){
  //   myservo1.write(i + 23); 
  //   myservo2.write(i + 59); 
  //   delay(1500);     
  //   }            // sets the servo position according to the scaled value

  delay(100);                           // waits for the servo to get there
}