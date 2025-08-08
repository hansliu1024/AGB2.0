You can run our tool on Windows using the command line interface. Our script runs with Python 3 and the NumPy library.

The input format for our tool comprises two types: one type tailored for regular LPN over the binary field, and another type for regular LPN over a larger field with field size |F|>2.

Assume script.py (our tool) is located in the C:\ directory.

==================Parameters==========
The order of command line parameters (n, k, t) is fixed, and we will explain their meanings.

n: number of samples; 
k: dimension of LPN; 
t: Hamming weight of a noise vector. 

======================================================================================================================== There is an input format for the regular LPN problem over the binary field. ======================

The input format for estimating the bit security of the regular LPN problem over the binary field is ``C:\AGBscript.py n=1024 k=652 t=57''

======================================================================================================================== There is an input format for the regular LPN problem over a larger field with field size |F|>2. ======================

The input format for estimating the bit security of the regular LPN problem over a larger field is ``C:\AGBscript.py n=1024 k=652 t=57 q''
