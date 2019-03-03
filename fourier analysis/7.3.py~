#!/usr/bin/python
import numpy as np
import pylab

def part_a():
    piano = np.loadtxt("piano.txt")
    trumpet = np.loadtxt("trumpet.txt")
    fft_piano = np.fft.rfft(piano)
    fft_trumpet = np.fft.rfft(trumpet)
    pylab.plot(fft_piano[:10000])
    pylab.plot(fft_trumpet[:10000])
    pylab.show()

def part_b():
    pass

def main():
    part_a()
    part_b()

if __name__ == "__main__":
    main()
