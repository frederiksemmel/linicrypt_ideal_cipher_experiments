# Modeling the ideal cipher in Linicrypt experiments

### Compression functions

This code lists all the compression functions from [PGV] and checks
for a collision structure.

The results are formatted in a grid, for convenience:
```
M=001        M=001        M=001        M=001
k=000        k=000        k=000        k=000
x=000        x=010        x=100        x=110
y=001        y=001        y=001        y=001
Degenerate   Degenerate   Degenerate   Degenerate

M=001        M=001        M=001        M=001
k=010        k=010        k=010        k=010
x=000        x=010        x=100        x=110
y=001        y=001        y=001        y=001
Degenerate   Degenerate   B            B

M=001        M=001        M=001        M=001
k=100        k=100        k=100        k=100
x=000        x=010        x=100        x=110
y=001        y=001        y=001        y=001
Degenerate   B            Degenerate   B

M=001        M=001        M=001        M=001
k=110        k=110        k=110        k=110
x=000        x=010        x=100        x=110
y=001        y=001        y=001        y=001
Degenerate   B            B            Degenerate

M=011        M=011        M=011        M=011
k=000        k=000        k=000        k=000
x=000        x=010        x=100        x=110
y=001        y=001        y=001        y=001
Degenerate   Degenerate   Degenerate   Degenerate

M=011        M=011        M=011        M=011
k=010        k=010        k=010        k=010
x=000        x=010        x=100        x=110
y=001        y=001        y=001        y=001
Degenerate   Degenerate   B            B

M=011        M=011        M=011        M=011
k=100        k=100        k=100        k=100
x=000        x=010        x=100        x=110
y=001        y=001        y=001        y=001
A            Secure       A            Secure

M=011        M=011        M=011        M=011
k=110        k=110        k=110        k=110
x=000        x=010        x=100        x=110
y=001        y=001        y=001        y=001
A            Secure       Secure       A

M=101        M=101        M=101        M=101
k=000        k=000        k=000        k=000
x=000        x=010        x=100        x=110
y=001        y=001        y=001        y=001
Degenerate   Degenerate   Degenerate   Degenerate

M=101        M=101        M=101        M=101
k=010        k=010        k=010        k=010
x=000        x=010        x=100        x=110
y=001        y=001        y=001        y=001
A            A            Secure       Secure

M=101        M=101        M=101        M=101
k=100        k=100        k=100        k=100
x=000        x=010        x=100        x=110
y=001        y=001        y=001        y=001
Degenerate   B            Degenerate   B

M=101        M=101        M=101        M=101
k=110        k=110        k=110        k=110
x=000        x=010        x=100        x=110
y=001        y=001        y=001        y=001
A            Secure       Secure       A

M=111        M=111        M=111        M=111
k=000        k=000        k=000        k=000
x=000        x=010        x=100        x=110
y=001        y=001        y=001        y=001
Degenerate   Degenerate   Degenerate   Degenerate

M=111        M=111        M=111        M=111
k=010        k=010        k=010        k=010
x=000        x=010        x=100        x=110
y=001        y=001        y=001        y=001
A            A            Secure       Secure

M=111        M=111        M=111        M=111
k=100        k=100        k=100        k=100
x=000        x=010        x=100        x=110
y=001        y=001        y=001        y=001
A            Secure       A            Secure

M=111        M=111        M=111        M=111
k=110        k=110        k=110        k=110
x=000        x=010        x=100        x=110
y=001        y=001        y=001        y=001
Degenerate   B            B            Degenerate
```
