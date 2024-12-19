#include <Servo.h>
#include <Stepper.h>
#define SERVO 43
#define POTENZIOMETRO A8

Servo myservo;  // create servo object to control a servo
// twelve servo objects can be created on most boards

int stepsPerRevolution = 2048;
int rpm = 10;
Stepper mystepper (stepsPerRevolution, 51, 49, 47, 45);

void setup() {
// Servo motore con potenziometro
  // myservo.attach(SERVO);  // attaches the servo on pin 9 to the servo object

// Motore stepper 
  mystepper.setSpeed(rpm);
}

void loop() {
// Servo motore con potenziometro
  // int input = analogRead(POTENZIOMETRO);            // reads the value of the potentiometer (value between 0 and 1023)
  // int angolo = map(input, 0, 1023, 0, 180);     // scale it to use it with the servo (value between 0 and 180)
  // myservo.write(angolo);                  // sets the servo position according to the scaled value
  // delay(15);                           // waits for the servo to get there

// Motore stepper 
  // mystepper.step(stepsPerRevolution);
  // delay(500);

  // mystepper.step(-stepsPerRevolution);
  // delay(500);

  int input = analogRead(POTENZIOMETRO);            
  int angolo = map(input, 0, 1023, 0, 180); 
  mystepper.step(angolo);
  delay(1000);

}