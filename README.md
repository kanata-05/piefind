# piefind
Find any sequence of digits in pi, uses the chudnovsky algorithm to calcuate pi (where's the fun in a pre-calcuated one)
Usage: ./pifind [-t <int>] -s <int>
-t value is the number of seconds until timeout (program stop) in order to avoid strain on the system, if left unspecified,
the program will run until the string is found.

-s value is the sequence of digits to find.
