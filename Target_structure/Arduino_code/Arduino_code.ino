#include <LiquidCrystal.h>
#include <Servo.h>
#define SERVO1 23 

// initialize the library by associating any needed LCD interface pin
// with the arduino pin number it is connected to
const int rs = 51, en = 47, d4 = 43, d5 = 41, d6 = 39, d7 = 37;
LiquidCrystal lcd(rs, en, d4, d5, d6, d7);

Servo servo;  // create servo object to control a servo
int offset = 62; // +- 3

int angolo = 40;

void setup() {
  // set up the LCD's number of columns and rows:
  lcd.begin(16, 2);
  // Print a message to the LCD.

  // Scrivi il primo messaggio sulla prima riga
  lcd.setCursor(0, 0); // Cursore a colonna 0, riga 0
  lcd.print("COMPTON");

  // Scrivi il secondo messaggio sulla seconda riga
  lcd.setCursor(0, 1); // Cursore a colonna 0, riga 1
  lcd.print("DESPERADOS ;)");
  
  servo.attach(SERVO1);  
}

void loop() {    

  servo.write(angolo + offset);      

  lcd.clear();

  lcd.setCursor(0, 0); // Cursore a colonna 0, riga 0
  lcd.print("Target: ");
  lcd.print(90);

  lcd.setCursor(0, 1); // Cursore a colonna 0, riga 0
  lcd.print("Detector: ");
  lcd.print(angolo);

  delay(500);

}