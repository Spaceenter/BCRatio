{
  Element::ReadMassTableFile("data/mass_table.dat");
  Element::ReadRadiusTableFile("data/radius_table.dat");
  
  Element* H  = new Element("Hydrogen",   1,  1, 63.04,  52.0, 42.8);
  Element* He = new Element("Helium",     4,  2,   0.0,   0.0,  0.0);
  Element* Li = new Element("Lithium",    7,  3,   0.0,   0.0,  0.0);
  Element* B  = new Element("Boron",      11, 5,   0.0,   0.0,  0.0);
  Element* C  = new Element("Carbon",    12,  6, 42.70,  85.8, 59.2); 
  Element* N  = new Element("Nitrogen",  14,  7, 37.99,  89.7, 61.1); 
  Element* O  = new Element("Oxygen",    16,  8, 34.24,  90.2, 61.3); 
  Element* Al = new Element("Aluminium", 26, 13, 24.01, 107.2, 69.7); 
  Element* Si = new Element("Silicon",   28, 14, 21.82, 108.4, 70.2); 
  Element* Ar = new Element("Argon",     39, 18, 19.55, 119.7, 75.7);
  Element* Xe = new Element("Xenon",    131, 54,  8.48, 172.1,100.8); 
  Element* Na = new Element("Sodium",    22, 11, 27.74, 102.6, 67.4); 
  Element* F  = new Element("Fluorine",  19,  9, 32.93,  94.4, 65.0); 
  Element* Fe = new Element("Iron",      56, 26,   0.0,   0.0,  0.0);
  Element* Pb = new Element("Lead",     207, 82,  6.37, 199.6,114.1);
}

