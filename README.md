# Modeling the ideal cipher in Linicrypt experiments

### Compression functions

This code lists all the compression functions from [PGV] and checks
for a collision structure.

The results are formatted in a grid, for convenience:
```
 M=001    M=001    M=001    M=001
0k=000   0k=000   0k=000   0k=000
0x=000   0x=010   0x=100   0x=110
0y=001   0y=001   0y=001   0y=001
 0,0,F    0,0,F    0,0,F    0,0,F
 0,0,B    0,0,B    0,0,B    0,0,B
                           
 M=001    M=001    M=001    M=001
0k=010   0k=010   0k=010   0k=010
0x=000   0x=010   0x=100   0x=110
0y=001   0y=001   0y=001   0y=001
 0,0,F    0,0,F    0,0,F    0,0,F
 0,0,B    0,0,B   Y0,0,B   Y0,0,B
                           
 M=001    M=001    M=001    M=001
0k=100   0k=100   0k=100   0k=100
0x=000   0x=010   0x=100   0x=110
0y=001   0y=001   0y=001   0y=001
 0,0,F    0,0,F    0,0,F    0,0,F
 0,0,B   Y0,0,B    0,0,B   Y0,0,B
                           
 M=001    M=001    M=001    M=001
0k=110   0k=110   0k=110   0k=110
0x=000   0x=010   0x=100   0x=110
0y=001   0y=001   0y=001   0y=001
 0,0,F    0,0,F    0,0,F    0,0,F
 0,0,B   Y0,0,B   Y0,0,B    0,0,B
                           
 M=011    M=011    M=011    M=011
0k=000   0k=000   0k=000   0k=000
0x=000   0x=010   0x=100   0x=110
0y=001   0y=001   0y=001   0y=001
 0,0,F    0,0,F   Y0,0,F   Y0,0,F
 0,0,B    0,0,B   Y0,0,B   Y0,0,B
                           
 M=011    M=011    M=011    M=011
0k=010   0k=010   0k=010   0k=010
0x=000   0x=010   0x=100   0x=110
0y=001   0y=001   0y=001   0y=001
 0,0,F    0,0,F    0,0,F    0,0,F
 0,0,B    0,0,B   Y0,0,B   Y0,0,B
                           
 M=011    M=011    M=011    M=011
0k=100   0k=100   0k=100   0k=100
0x=000   0x=010   0x=100   0x=110
0y=001   0y=001   0y=001   0y=001
Y0,0,F    0,0,F   Y0,0,F    0,0,F
 0,0,B    0,0,B    0,0,B    0,0,B
                           
 M=011    M=011    M=011    M=011
0k=110   0k=110   0k=110   0k=110
0x=000   0x=010   0x=100   0x=110
0y=001   0y=001   0y=001   0y=001
Y0,0,F    0,0,F    0,0,F   Y0,0,F
 0,0,B    0,0,B    0,0,B    0,0,B
                           
 M=101    M=101    M=101    M=101
0k=000   0k=000   0k=000   0k=000
0x=000   0x=010   0x=100   0x=110
0y=001   0y=001   0y=001   0y=001
 0,0,F   Y0,0,F    0,0,F   Y0,0,F
 0,0,B   Y0,0,B    0,0,B   Y0,0,B
                           
 M=101    M=101    M=101    M=101
0k=010   0k=010   0k=010   0k=010
0x=000   0x=010   0x=100   0x=110
0y=001   0y=001   0y=001   0y=001
Y0,0,F   Y0,0,F    0,0,F    0,0,F
 0,0,B    0,0,B    0,0,B    0,0,B
                           
 M=101    M=101    M=101    M=101
0k=100   0k=100   0k=100   0k=100
0x=000   0x=010   0x=100   0x=110
0y=001   0y=001   0y=001   0y=001
 0,0,F    0,0,F    0,0,F    0,0,F
 0,0,B   Y0,0,B    0,0,B   Y0,0,B
                           
 M=101    M=101    M=101    M=101
0k=110   0k=110   0k=110   0k=110
0x=000   0x=010   0x=100   0x=110
0y=001   0y=001   0y=001   0y=001
Y0,0,F    0,0,F    0,0,F   Y0,0,F
 0,0,B    0,0,B    0,0,B    0,0,B
                           
 M=111    M=111    M=111    M=111
0k=000   0k=000   0k=000   0k=000
0x=000   0x=010   0x=100   0x=110
0y=001   0y=001   0y=001   0y=001
 0,0,F   Y0,0,F   Y0,0,F    0,0,F
 0,0,B   Y0,0,B   Y0,0,B    0,0,B
                           
 M=111    M=111    M=111    M=111
0k=010   0k=010   0k=010   0k=010
0x=000   0x=010   0x=100   0x=110
0y=001   0y=001   0y=001   0y=001
Y0,0,F   Y0,0,F    0,0,F    0,0,F
 0,0,B    0,0,B    0,0,B    0,0,B
                           
 M=111    M=111    M=111    M=111
0k=100   0k=100   0k=100   0k=100
0x=000   0x=010   0x=100   0x=110
0y=001   0y=001   0y=001   0y=001
Y0,0,F    0,0,F   Y0,0,F    0,0,F
 0,0,B    0,0,B    0,0,B    0,0,B
                           
 M=111    M=111    M=111    M=111
0k=110   0k=110   0k=110   0k=110
0x=000   0x=010   0x=100   0x=110
0y=001   0y=001   0y=001   0y=001
 0,0,F    0,0,F    0,0,F    0,0,F
 0,0,B   Y0,0,B   Y0,0,B    0,0,B
```

The last two rows indicate whether the program has a collision structure of that type.
In the case of a single query program, there is only type forward and type backward.

There are 18 programs of type F and 18 of type B.
12 are secure, meaning they have no collision structure.

### Finding interesting examples for 3 inputs and 2 queries

Below is a small piece of the output if we search for programs that have <= 2 collision structures.
There are 12 possible collision structures for programs of this type.
```
 M=00001    M=00001    M=00001    M=00001    M=00001    M=00001    M=00001    M=00001
0k=00100   0k=00100   0k=00100   0k=00100   0k=00100   0k=00100   0k=00100   0k=00100
0x=01000   0x=01000   0x=01000   0x=01000   0x=01000   0x=01000   0x=01000   0x=01000
0y=00010   0y=00010   0y=00010   0y=00010   0y=00010   0y=00010   0y=00010   0y=00010
1k=01010   1k=01010   1k=01010   1k=01010   1k=01010   1k=01010   1k=01010   1k=01010
1x=10000   1x=10010   1x=10100   1x=10110   1x=11000   1x=11010   1x=11100   1x=11110
1y=00001   1y=00001   1y=00001   1y=00001   1y=00001   1y=00001   1y=00001   1y=00001
 01,0,FF    01,0,FF    01,0,FF    01,0,FF    01,0,FF    01,0,FF    01,0,FF    01,0,FF
Y01,0,FB   Y01,0,FB   Y01,0,FB   Y01,0,FB   Y01,0,FB   Y01,0,FB   Y01,0,FB   Y01,0,FB
 01,0,BF    01,0,BF    01,0,BF    01,0,BF    01,0,BF    01,0,BF    01,0,BF    01,0,BF
Y01,0,BB   Y01,0,BB   Y01,0,BB   Y01,0,BB   Y01,0,BB   Y01,0,BB   Y01,0,BB   Y01,0,BB
 10,0,FF    10,0,FF    10,0,FF    10,0,FF    10,0,FF    10,0,FF    10,0,FF    10,0,FF
 10,0,FB    10,0,FB    10,0,FB    10,0,FB    10,0,FB    10,0,FB    10,0,FB    10,0,FB
 10,0,BF    10,0,BF    10,0,BF    10,0,BF    10,0,BF    10,0,BF    10,0,BF    10,0,BF
 10,0,BB    10,0,BB    10,0,BB    10,0,BB    10,0,BB    10,0,BB    10,0,BB    10,0,BB
 01,1,F     01,1,F     01,1,F     01,1,F     01,1,F     01,1,F     01,1,F     01,1,F
 01,1,B     01,1,B     01,1,B     01,1,B     01,1,B     01,1,B     01,1,B     01,1,B
 10,1,F     10,1,F     10,1,F     10,1,F     10,1,F     10,1,F     10,1,F     10,1,F
 10,1,B     10,1,B     10,1,B     10,1,B     10,1,B     10,1,B     10,1,B     10,1,B

...

 M=00011    M=00011    M=00011    M=00011    M=00011    M=00011    M=00011    M=00011
0k=01100   0k=01100   0k=01100   0k=01100   0k=01100   0k=01100   0k=01100   0k=01100
0x=00100   0x=00100   0x=00100   0x=00100   0x=00100   0x=00100   0x=00100   0x=00100
0y=00010   0y=00010   0y=00010   0y=00010   0y=00010   0y=00010   0y=00010   0y=00010
1k=11000   1k=11000   1k=11000   1k=11000   1k=11000   1k=11000   1k=11000   1k=11000
1x=00000   1x=00010   1x=00100   1x=00110   1x=01000   1x=01010   1x=01100   1x=01110
1y=00001   1y=00001   1y=00001   1y=00001   1y=00001   1y=00001   1y=00001   1y=00001
 01,0,FF    01,0,FF    01,0,FF    01,0,FF    01,0,FF    01,0,FF    01,0,FF    01,0,FF
 01,0,FB    01,0,FB    01,0,FB    01,0,FB    01,0,FB    01,0,FB    01,0,FB    01,0,FB
 01,0,BF    01,0,BF    01,0,BF    01,0,BF    01,0,BF    01,0,BF    01,0,BF    01,0,BF
 01,0,BB    01,0,BB    01,0,BB    01,0,BB    01,0,BB    01,0,BB    01,0,BB    01,0,BB
 10,0,FF    10,0,FF    10,0,FF    10,0,FF    10,0,FF    10,0,FF    10,0,FF    10,0,FF
Y10,0,FB    10,0,FB    10,0,FB    10,0,FB    10,0,FB    10,0,FB   Y10,0,FB   Y10,0,FB
 10,0,BF    10,0,BF    10,0,BF    10,0,BF    10,0,BF    10,0,BF    10,0,BF    10,0,BF
 10,0,BB    10,0,BB    10,0,BB    10,0,BB    10,0,BB    10,0,BB   Y10,0,BB   Y10,0,BB
 01,1,F     01,1,F     01,1,F     01,1,F     01,1,F     01,1,F     01,1,F     01,1,F
 01,1,B     01,1,B     01,1,B     01,1,B     01,1,B     01,1,B     01,1,B     01,1,B
 10,1,F     10,1,F     10,1,F     10,1,F     10,1,F     10,1,F     10,1,F     10,1,F
Y10,1,B    Y10,1,B     10,1,B     10,1,B     10,1,B     10,1,B     10,1,B     10,1,B

...
```

These combinations of collision structures occur.
```
000000000000: 76668
000000000001: 3696
000000000010: 5568
000000000011: 4524
000000010001: 2760
000000100010: 2544
000000110011: 1176
000001000001: 2760
000001010000: 11208
000001010001: 684
000010000010: 2544
000010100000: 12360
000010100010: 636
000010101100: 1392
000011000011: 1176
000011110011: 240
010000000100: 5520
010000100110: 2784
010010100100: 2352
010100000000: 11268
010100000011: 348
010100000100: 1368
010100010001: 1392
010100100010: 1272
010100110011: 696
010100110111: 1392
010101010000: 2088
010101010001: 696
010101010100: 1392
010110100000: 1848
010110100010: 636
010110100100: 1176
010111110011: 348
100000001000: 5520
100010001010: 2784
100010101000: 2352
101000000000: 11268
101000000011: 348
101000001000: 1368
101001000001: 1392
101001010000: 2088
101001010001: 696
101001011000: 1392
101010000010: 1272
101010100000: 1848
101010100010: 636
101010101000: 1176
101011000011: 696
101011001011: 1392
101011110011: 348
110000001100: 2784
110010101100: 1392
110010101110: 1392
111100001100: 696
111101011100: 696
111110101100: 696
111111111111: 696
```

Note, that no program has only a single collision structure with 2 different queries.

### Better compression ratio: 3 blocks using 2 queries

As an experiment, I looked for such a linicrypt program without a collision structure.
This is the first such a program, sorting the programs lexicographically.

```
 M=000101
0k=001000
0x=010000
0y=000010
1k=100000
1x=010110
1y=000001
```
