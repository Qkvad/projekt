# Diskontinuirana metoda konacnih elemenata za elipticku jednadzbu 

U ovom radu je prezentirana diskontinuirana Galerkinova metoda konacnih elemenata za tri razlicita primjera eliptickih jednadzbi koja su razlicite glatkoce i na drugacijim mrezama. 

---

Kod se moze naci u ``src`` direktoriju, a seminarski rad koji opisuje teoriju i kod se moze naci u direktoriju ``seminar``.

---

U ``src/projekt.cc`` treba definirati problem koji zelimo rijesiti i parametre za njega, te na jasno naznacenom mjestu u kodu treba staviti odgovarajuci put do mreze (ovdje je to onako kako to izgleda na virtualnoj masini).

Parametre za primjer 1 se moze naci u ``src/parameterA.hh`` te odgovarajuci problem je dan u ``src/problemA.hh``.
Metodu se moze pokrenuti uz npr sljedece argumente komandne linije 
```2 simplex 1 NIPG 3 1```
Parametre za primjer 2 se moze naci u ``src/parameterB.hh`` te odgovarajuci problem je dan u ``src/problemB.hh``.
Metodu se moze pokrenuti uz npr sljedece argumente komandne linije 
```2 simplex 1 NIPG 3 1```
