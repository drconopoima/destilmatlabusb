# destilmatlabusb
Versatile Matlab thermodynamics package that resolves Distillation columns by using shortcut methods and rigorous methods, enabled with a Matlab GUI, and powered by fully adaptable/programmable equations of state (EdE) that can be overwritten with own parameters or which methods can be replaced with in-built Matlab functions to model different equations. Comes pre-packaged with Peng-Robinson, Soave-Redlich-Kwong and Van der Waals EoS

## Getting Started

Get the full DESTILMATLAB folder into your Matlab working directory and change
to the internal of the directory.

The project was coded using the Object-Oriented programming MATLAB syntax.

Initially, you will want to retrieve from a database the properties of pure substances

The Class object that is used to store properties is @SUSTANCIA. The method that
queries the DIPPR database is get_properties('Substance Name').

You can create a substance by passing the string with the IUPAC name to the class
Sustancia.

```
$ Sustancia('ethane')
```

The resolvers for currents and for distillation towers typically accept
conditions in the current at the feed. For which you need to define the whole
composition of the feed. You can do so in a List:

```
sust = [Sustancia('ethane'), Sustancia('Propane'), Sustancia('Butane'), Sustancia('Pentane'),
Sustancia('Hexane')];
```

The composition on a feed can also be passed as a list. Lets assume there are
several feeds and define both with concentrations:

```
conc1 = [0.5342, 0.3163, 0.0612, 0.0701, 0.0183];
```
```
conc2 = [0.2236, 0.5194, 0.1119, 0.1230, 0.0221];
```

If the interaction parameters kij are taken as ideal:

```
kij = 0;
```

Then to fully define the feeds, we create the following mixtures

```
mix1 = Mezcla(sust, conc1, kij);
mix2 = Mezcla(sust, conc2, kij);
```

If we wanted to use the Equation of State formulated by Peng and Robinson:

```
EdE = PREdE();
```

If we wanted to use mixture rules formulated based on Van der Waals:



Please check the resolved examples in each folder.
