import numpy as np

def bit_rate(Mpix, FR, bit_depth):
    # Mpix: megapixels
    # FR: frame rate in Hz
    # bit_depth: in bits
    bit_rate = Mpix*1e6*FR*bit_depth
    return bit_rate

def total_data(bit_rate, time):
    # bit_rate: in bits per second
    # time: in seconds
    total_data = bit_rate*time
    return total_data

print('# Data output calculator')
Mpix = float(input('Camera megapixels (eg 67): ').strip() or '67')
FR = float(input('Frame rate in Hz (eg 90): ').strip() or '90')
bit_depth = float(input('Bit depth (eg 10): ').strip() or '10')
time = float(input('Time in min (eg 60): ').strip() or '60')

bit_rate = bit_rate(Mpix, FR, bit_depth)
total_data = total_data(bit_rate, time*60)

print('\n--------------')
print('Bit rate = ' + str(bit_rate/8e9) + ' Gb/s')
print('Total data = ' + str(total_data/8e12) + ' Tb in ' + str(time) + ' min')
print('--------------\n')

