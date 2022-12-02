# from stack overflow: ~/questions/21548464/how-to-write-a-makefile-to-compile-a-simple-c-program 
CC				= gcc
CC_FLAGS 	= -mavx -mfma -mavx2 -O3 -std=c99 
RM 				= rm -f

default: all

all: kernel

kernel: kernel_driver.x
	./kernel_driver.x

kernel_driver.x:
	$(CC) $(CC_FLAGS) -o kernel_driver.x kernel/kernel_driver.c kernel/kernel.c -lm


clean:
	rm -rf *.x *.S
