#include <LiquidCrystal.h>
#include <Servo.h>
#define SERVO1 23 
#define SERVO2 29

// initialize the library by associating any needed LCD interface pin
// with the arduino pin number it is connected to
const int rs = 51, en = 47, d4 = 43, d5 = 41, d6 = 39, d7 = 37;
LiquidCrystal lcd(rs, en, d4, d5, d6, d7);

Servo alto;  // create servo object to control a servo
Servo basso;  // create servo object to control a servo
int offset_alto = 19; //23 x0 -- 
int offset_basso = 55; //59 x0 --

int angolo_alto = 90;
int angolo_basso = 90;

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
  
  alto.attach(SERVO1);  
  basso.attach(SERVO2);  
}

void loop() {    

  alto.write(offset_alto + angolo_alto);     
  basso.write(offset_basso + angolo_basso);

  lcd.clear();

  lcd.setCursor(0, 0); // Cursore a colonna 0, riga 0
  lcd.print("Target: ");
  lcd.print(angolo_alto + 90);

  lcd.setCursor(0, 1); // Cursore a colonna 0, riga 0
  lcd.print("Detector: ");
  lcd.print(angolo_basso);

  delay(500); 
}