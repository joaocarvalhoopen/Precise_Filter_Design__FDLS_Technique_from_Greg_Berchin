###############################################################################
# Project name: An implementation of a Precise Filter Design.                 #
# Author: João Nuno Carvalho                                                  #
# Date:   2020.07.04                                                          #
# Method Author: Greg Berchin                                                 #
#                                                                             #
# Description: This is a simple implementation in Python of the FDLS          #
#              (Frequency Domain Least Squares) technique described in        #
#              chapter 7 of the book Streamlining Digital Signal Processing   #
# 2th Ed. . The same technique was also described in the Master thesis of the #
# author of the technique, Greg Berchin.                                      #
# The technique is a method of filter construction. By starting from the      #
# specification of the frequency response that we requirer (attenuation and   #
# phase-shift by frequency plots), we can design a FIR or a IIR filter that   #
# is in principal more precise that the normally designed FIR or IIR filter   #
# and it doesn't need to be a linear phase filter like the normal FIR designs.#
# This technique has many different applications, but the author gives the    #
# following example ... Imagine that you were an Electrical Engineer given    #
# the task of replacing an analog filter circuit or section with digital      #
# filter. With this technique you could start by testing and characterizing   #
# the analog circuit, obtain the frequency response of the analog circuit,    #
# it's attenuation with frequency and it's phase-shift with frequency and     #
# then using the technique, design a FIR or a IIR digital filter that would   #
# replace the original analog circuit.                                        #
#                                                                             #
# Note: See the book for more details.                                        #
#                                                                             #
# License: MIT Opensource License.                                            #
#                                                                             #
###############################################################################

#######
# Inputs:
#    list_frequencies
#    list_attenuations
#    list_phase_shift

#######
# Outputs:
#    list_FIR_filter_values

#######
# Signal input data:
#    test_signal_input

#######
# Signal output data:
#    test_signal_output

#######
# Equations and steps of the FIR and IIR filter design:


'''

###########
# Equations

 =>Filter in polynomial form:

     Y(Z)    b_0 + b_1*Z^-1 + ...... + b_n*Z^-N
     ---- = ------------------------------------            (1)
     U(Z)     1  + a_1*Z^-1 + ...... + a_d*Z^-D

 
 =>Matrix system that we want to solve:

    Y = X * THETA                                           (2)

 
 =>To fill the Y column vector of equation (2):
 
        ---              --- 
        | A_1 * cos(phi_1) |
    Y = | A_2 * cos(phi_2) |                                (3)
        |      ......      |
        | A_m * cos(phi_m) |
        ---              ---


=>To fill de X matrix of equation (2):

                                                                     ---   ---
 ---      ---   ---                                              --- |  a_1  |
 | y_1( 0 ) |   | -y_1( -1 ) ... -y_1( -D ) u_1( 0 ) ... u_1( -N ) | |  ...  |
 | y_2( 0 ) |   | -y_2( -1 ) ... -y_2( -D ) u_2( 0 ) ... u_2( -N ) | |  ...  |
 |   ...    | = |   ...            ...        ...         ...      | |  a_D  |    (4)
 |   ...    |   |   ...            ...        ...         ...      | |  b_0  | 
 | y_M( 0 ) |   | -y_M( -1 ) ... -y_M( -D ) u_M( 0 ) ... u_M( -N ) | |  ...  |
 ---      ---   ---                                              --- |  ...  |
                                                                     |  b_N  |
                                                                     ---   ---

 =>How to fill one line of matrix X, for the case that D = 2 and N = 2.

   Starting with:

     -y_1( -1 ) ... -y_1( -D ) u_1( 0 ) ... u_1( -N )                    (5)

   it becomes:

     -y_1( -1 ) -y_1( -2 ) u_1( 0 ) u_1( -1 ) u_1( -2 )                  (6)

   and that then becomes:

     -y_1( -1 ) = A_1 * cos((-1)*W_1*T_s + phi_1) = -A_1 * cos(-W_1*T_s + phi_1)

     -y_1( -2 ) = A_1 * cos((-2)*W_2*T_s + phi_2) = -A_1 * cos(-W_2*T_s + phi_2)

      u_1( 0 ) = cos((0)*W_1*T_s) = cos(0*w_1*T_s) = 1                           (7 as an all)

      u_1( -1 ) = cos((-1)*W_1*T_s) = cos(-W_1*T_s)

      u_1( -2 ) = cos((-2)*W_1*T_s) = cos(-2*W_1*T_s)


 =>don't forget that:

      w_M = 2*PI*freq_M                                   (8) 

   in witch the M is the max number of data points but in here
   corresponds to one data point of the M possible, 
   that is the frequency of one data point in the CSV table.


 =>And have in mind that T_s is the period of
   the sample_rate of the signal, that is:

   for a sample_rate = 1000 samples / second  it will be
    
     T_s = 1 / sample_rate                                (9)
   
   that is T_s = 0.001 = 10^-3 seconds 


 =>Pseudo-Inverse calculation to obtain the column vector THETA of equation (2),
   this vector will have the two types of coefficients for the filter:

    THETA ≈ Inverse(Transpose(X) * X) * Transpose(X) * Y        (10)


'''

# Steps of the FIR filter design:
#
# 0. Read configuration file with filter program parameters. It's a CSV.
#
# 1. Start with a file containing the frequency response that we would like
#    that our filter to have. This CSV file should have many data points. 
#    The CSV file has the following structure: <br>
#
#    frequency, amplitude (in linear scale), phase (in radians) <br>
#
#    An interpolation function could be optionally added to generate more
#    intermediate points.
#    In here, we will have to choose the M, the number of data points (tuple)
#    that we will be using. The M will dictate the number of rows of the matrix X.
#
# 2. Choose the N=? and D=? values. The N is the degree of the polynomial in
#    the numerator and the D is the degree of the polynomial in the
#    denominator of equation (1). Only one N or D can be zero at the same time,
#    but they can be a zero or a positive number of any value. The greater the N
#    and D number, the larger the filter will be that we are constructing and the
#    more demanding in terms of computational cost it will be to run it.
#
# 3. Choose the sample_rate of the digital signal to choose the T_s, that is the
#    period of the sample rate that is calculated with equation (9).
#
# 4. Fill the y_m sequence of cosines for the length of D.
#
# 5. Fill the u_m sequence of cosines for the length of N+1 (Note: In here is
#    plus one terms).
#
# 6. Join the columns of y_m concatenated with u_m to make the matrix X of
#    matrix equation X (2).
#
# 7. Create the column vector Y of equation (3).
#
# 8. Calculate the Pseudo-Inverse the vector THETA. It will have the filter
#    coefficients. Depending on the number chosen for N and D, a different
#    type of filter topology will be made.
#
# 9. Construct the filter.
#
# 10. Test the filter frequency response with a chirp signal between the
#     frequency limits. 
#
#


#######
# Steps of the FIR filter test:
#
# 1- For the each of different input test signals generate them.
# 2- Plot the spectrogram of input signal, magnitude and phase-shift using the FFT.
# 3- Apply the designed filter to the input signal and obtain the output signal.
# 4- Plot the spectrogram of output signal, magnitude and phase-shift using the FFT. 
#
# List of test signals to be generated:
#
#     1. Simple sine or cosine functions with and without phase shifts.
#     2. Two crafted sinusoids combined.
#     3. N crafted sinusoids combined.
#     4. Chirp from freq A to freq B can also be called a frequency sweep.
#     5. One Dirac impulse.
#     6. One step cycle function.    
#     7. One square cycle function.
#     8. Square wave function.
#     9. Uniform noise from freq A to freq B.
#

#############################
#############################


import csv
import math


import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import inv
from numpy.linalg import pinv


# Inputs:
list_frequencies  = []
list_attenuations = []
list_phase_shift  = []

# Outputs:
list_impulse_responce  = []
list_FIR_filter_values = []


##########
# 1. Sine Signal

def generate_sine(sample_rate, signal_duration, freq, phase, amplitude):
    freq      = float(freq)
    phase     = float(phase)
    amplitude = float(amplitude)
    start = 0.0
    stop = signal_duration
    step = 1.0 / sample_rate
    t = np.arange(start, stop, step, dtype='float')
    s = amplitude*np.sin(phase + 2*np.pi*freq*t)
    len_s = s.size
    return (t, s, len_s)

def plot_sine(t, s, freq, phase):
    plt.plot(t, s)
    plt.xlabel('time (s)')
    plt.ylabel('amplitude')
    plt.title('Sine signal freq. ' + str(freq) + ' Hz' + ' phase ' + "{0:.2f}".format(phase) + ' Rad'  )
    plt.grid(True)
    #plt.savefig("test.png")
    plt.show()

def test_sine_signal():
    # Generate sine 1.
    sample_rate = 40000
    freq        = 500.0
    phase       = 0.0
    amplitude   = 1.0
    signal_duration = 0.01
    t, s, len_s = generate_sine(sample_rate, signal_duration, freq, phase, amplitude)
    plot_sine(t, s, freq, phase)

    # Generate sine 2.
    sample_rate = 40000
    freq        = 500.0
    phase       = np.pi / 2.0
    amplitude   = 1.0
    signal_duration = 0.01
    t, s, len_s = generate_sine(sample_rate, signal_duration, freq, phase, amplitude)
    plot_sine(t, s, freq, phase)
    plot_spectrum(t, s, sample_rate)

#########
# 2. Two crafted sinusoids combined.

def generate_2_sines(sample_rate, signal_duration, freq_1, phase_1, amplitude_1, freq_2, phase_2, amplitude_2):
    freq_1      = float(freq_1)
    phase_1     = float(phase_1)
    amplitude_1 = float(amplitude_1)
    freq_2      = float(freq_2)
    phase_2     = float(phase_2)
    amplitude_2 = float(amplitude_2)
    t_1, s_1, len_s_1 = generate_sine(sample_rate, signal_duration, freq_1, phase_1, amplitude_1)
    t_2, s_2, len_s_2 = generate_sine(sample_rate, signal_duration, freq_2, phase_2, amplitude_2)
    s_output = (s_1 + s_2) / 2.0
    len_s = s_output.size
    return (t_1, s_output, len_s)

def plot_2_sines(t, s, freq_1, phase_1, freq_2, phase_2 ):
    plt.plot(t, s)
    plt.xlabel('time (s)')
    plt.ylabel('amplitude')
    plt.title('Two sines signals freq_a ' + str(freq_1) + ' Hz' + ' phase_a ' + "{0:.2f}".format(phase_1) + ' Rad' +
              ' freq_b ' + str(freq_2) + ' Hz' + ' phase_b ' + "{0:.2f}".format(phase_2) + ' Rad' )
    plt.grid(True)
    #plt.savefig("test.png")
    plt.show()

def test_2_sines_signal():
    sample_rate = 40000
    freq_1      = 500.0
    phase_1     = 0.0
    amplitude_1 = 1.0
    freq_2      = 750.0
    phase_2     = np.pi*(2/3.0)
    amplitude_2 = 1.0
    signal_duration = 0.01
    t, s, len_s = generate_2_sines(sample_rate, signal_duration, freq_1, phase_1, amplitude_1,
                                   freq_2, phase_2, amplitude_2)
    plot_2_sines(t, s, freq_1, phase_1, freq_2, phase_2)
    plot_spectrum(t, s, sample_rate)

#########
# 3. N crafted sinusoids combined.

def generate_N_sines(sample_rate, signal_duration, freq_list, phase_list, amplitude_list):
    len_a = len(freq_list)
    len_b = len(phase_list)
    len_c = len(amplitude_list)
    if ((len_a != len_b != len_c) and len_a > 0 and len_b > 0 and len_c > 0):
        print("ERROR: In function generate_N_sines() parameters freq_list, phase_list and amplitude_list are not of equal size and greater then zero.")
        exit(-1)
    t_output = None
    accumulator = None
    len_seq = len(freq_list)
    for i in range(len_seq):
        freq      = float(freq_list[i])
        phase     = float(phase_list[i])
        amplitude = float(amplitude_list[i]) 
        t, s, len_s = generate_sine(sample_rate, signal_duration, freq, phase, amplitude)
        if t_output is None:
            t_output = t
        accumulator = s if (accumulator is None) else accumulator + s
    s_output = accumulator / float(len_seq) 
    len_s = s_output.size
    return (t_output, s_output, len_s)

def plot_N_sines(t, s, freq_list, phase_list):
    plt.plot(t, s)
    plt.xlabel('time (s)')
    plt.ylabel('amplitude')
    str_freq_phase = ["freq {0:.2f} Hz, phase {1:.2f} Rad".format(pair[0], pair[1]) for pair in zip(freq_list, phase_list)]    
    plt.title('N sines ' + "".join(str_freq_phase))
    plt.grid(True)
    #plt.savefig("test.png")
    plt.show()

def test_N_sines_signal():
    sample_rate     = 40000
    freq_list       = [500.0, 750.0, 950.0]
    phase_list      = [0.0, np.pi/2.0, np.pi*(2/3)]
    amplitude_list  = [1.0, 1.0, 1.0]
    signal_duration = 0.01
    t, s, len_s = generate_N_sines(sample_rate, signal_duration, freq_list, phase_list, amplitude_list)
    plot_N_sines(t, s, freq_list, phase_list)
    plot_spectrum(t, s, sample_rate)

##########
# 4. Chirp Signal

def generate_chirp(sample_rate, chirp_duration, start_freq, end_freq):
    f_0 = start_freq
    f_1 = end_freq
    start = 0.0
    stop = chirp_duration
    step = 1.0 / sample_rate
    t = np.arange(start, stop, step, dtype='float')
    phase = 0.0
    chirp_period = chirp_duration # 1 / 100.0 #1.0
    k = (f_1 - f_0) / chirp_period
    s = np.sin(phase + 2*np.pi * ( f_0*t + (k/2)*np.square(t)) )
    len_s = s.size
    return (t, s, len_s)

def plot_chirp(t, s, start_freq, end_freq):
    plt.plot(t, s)
    plt.xlabel('time (s)')
    plt.ylabel('amplitude')
    plt.title('Chirp  from ' + str(start_freq) + 'Hz  to ' + str(end_freq) + 'Hz' )
    plt.grid(True)
    #plt.savefig("test.png")
    plt.show()

def test_chirp_signal():
    sample_rate = 40000
    start_freq = 500
    end_freq   = 5000
    # chirp_duration = 1.00
    chirp_duration = 0.01
    t, s, len_s = generate_chirp(sample_rate, chirp_duration, start_freq, end_freq)
    plot_chirp(t, s, start_freq, end_freq)
    plot_spectrum(t, s, sample_rate)

##########
# 5. One Dirac impulse.

def generate_dirac_impulse(sample_rate, buffer_duration, dirac_impulse_offset):
    start = 0.0
    stop = buffer_duration
    step = 1.0 / sample_rate
    t = np.arange(start, stop, step, dtype='float')
    s = np.zeros(len(t), dtype='float')
    s[dirac_impulse_offset] = 1.0
    len_s = s.size
    return (t, s, len_s)

def plot_dirac_impulse(t, s, dirac_impulse_offset):
    plt.plot(t, s)
    plt.xlabel('time (s)')
    plt.ylabel('amplitude')
    plt.title('Dirac impulse offset ' + str(dirac_impulse_offset) + ' samples' )
    plt.grid(True)
    #plt.savefig("test.png")
    plt.show()

def test_dirac_impulse_signal():
    sample_rate = 40000
    # chirp_duration = 1.00
    buffer_duration = 0.01
    dirac_impulse_offset = 20  # Samples
    t, s, len_s = generate_dirac_impulse(sample_rate, buffer_duration, dirac_impulse_offset)
    plot_dirac_impulse(t, s, dirac_impulse_offset)
    plot_spectrum(t, s, sample_rate)

###########
# 6. One step function.

def generate_step(sample_rate, buffer_duration, step_impulse_offset):
    start = 0.0
    stop  = buffer_duration
    step  = 1.0 / sample_rate
    t = np.arange(start, stop, step, dtype='float')
    s = np.zeros(len(t), dtype='float')
    for i in range(0, len(s)):
        if i >= step_impulse_offset:
            s[i] = 1
    len_s = s.size
    return (t, s, len_s)

def plot_step(t, s, step_impulse_offset):
    plt.plot(t, s)
    plt.xlabel('time (s)')
    plt.ylabel('amplitude')
    plt.title('One step impulse offset ' + str(step_impulse_offset) + ' samples' )
    plt.grid(True)
    #plt.savefig("test.png")
    plt.show()

def test_step_signal():
    sample_rate = 40000
    # chirp_duration = 1.00
    buffer_duration = 0.01
    step_impulse_offset = 20  # Samples
    t, s, len_s = generate_step(sample_rate, buffer_duration, step_impulse_offset)
    plot_step(t, s, step_impulse_offset)
    plot_spectrum(t, s, sample_rate)

###########
# 7. One square function.

def generate_square(sample_rate, buffer_duration, square_impulse_offset):
    start = 0.0
    stop  = buffer_duration
    step  = 1.0 / sample_rate
    t = np.arange(start, stop, step, dtype='float')
    s = np.zeros(len(t), dtype='float')
    for i in range(0, len(s)):
        if (i >= square_impulse_offset and i <= 2*square_impulse_offset ):
            s[i] = 1
    len_s = s.size
    return (t, s, len_s)

def plot_square(t, s, square_impulse_offset):
    plt.plot(t, s)
    plt.xlabel('time (s)')
    plt.ylabel('amplitude')
    plt.title('Square impulse offset_1 ' + str(square_impulse_offset) + ', samples offset_2 ' +
              str(2*square_impulse_offset+1) + ' samples' )
    plt.grid(True)
    #plt.savefig("test.png")
    plt.show()

def test_square_signal():
    sample_rate           = 40000
    buffer_duration       = 0.01
    square_impulse_offset = 20  # Samples
    t, s, len_s = generate_square(sample_rate, buffer_duration, square_impulse_offset)
    plot_square(t, s, square_impulse_offset)
    plot_spectrum(t, s, sample_rate)

#########
# 8. Square wave function.

def generate_square_wave(sample_rate, buffer_duration, square_wave_width):
    start = 0.0
    stop  = buffer_duration
    step  = 1.0 / sample_rate
    t = np.arange(start, stop, step, dtype='float')
    s = np.zeros(len(t), dtype='float')
    counter = square_wave_width
    flag_val = False 
    for i in range(0, len(s)):
        if counter == 0:
            flag_val = not flag_val
            counter = square_wave_width
        s[i] = 0.0 if (flag_val == False) else 1.0
        counter -= 1
    len_s = s.size
    return (t, s, len_s)

def plot_square_wave(t, s, square_wave_width):
    plt.plot(t, s)
    plt.xlabel('time (s)')
    plt.ylabel('amplitude')
    plt.title('Square wave width ' + str(square_wave_width) + ' samples' )
    plt.grid(True)
    #plt.savefig("test.png")
    plt.show()

def test_square_wave_signal():
    sample_rate      = 40000
    buffer_duration  = 0.01   # 1
    square_wave_width = 20    # Samples
    t, s, len_s = generate_square_wave(sample_rate, buffer_duration, square_wave_width)
    plot_square_wave(t, s, square_wave_width)
    plot_spectrum(t, s, sample_rate)

#########
# 9. Uniform noise from freq A to freq B.

def generate_noise(sample_rate, signal_duration, start_freq, end_freq):
    # We use the max sample_rate for the fft sample size, so that all
    # frequencies up until the max sample rate can be used.
    # The generated signal duration is always one second, with 
    # the buffer size of the sample rate.

    # Note: Inspired on the post of stackoverflow:
    #       How to generate noise in frequency range with numpy?
    # https://stackoverflow.com/questions/33933842/how-to-generate-noise-in-frequency-range-with-numpy

    # In here we determine the frequencies for sample_rate or fft buffer_size.
    fft_freq = np.abs(np.fft.fftfreq(sample_rate, 1 / sample_rate))
 
    chosen_freq_bins = np.zeros(sample_rate)
    for i in range(0, len(fft_freq)):
        if (fft_freq[i] >= start_freq and fft_freq[i] <= end_freq):
            chosen_freq_bins[i] = int(1)
    ifft_buffer = np.array(chosen_freq_bins, dtype=complex)
    N_positive = (len(ifft_buffer) - 1) // 2
    # Note: The amplitude is always 1 (one), the phase is what can change
    #       randomly.
    # Phase in the interval [0, 2*Pi].
    phase_value = np.random.rand(N_positive) * 2 * np.pi
    complex_value = np.cos(phase_value) + 1j * np.sin(phase_value)
    # The first half of the buffer has the phase.
    ifft_buffer[ 1: N_positive + 1] *= complex_value
    # The second half of the buffer has complex conjugate of the phase. 
    ifft_buffer[ -1 : -1-N_positive : -1] = np.conj(ifft_buffer[ 1: N_positive + 1])
    # Get the real part of the iFFT of the input buffer for the generated
    # noise for the selected frequencies.
    output_buffer = np.fft.ifft(ifft_buffer).real
    
    if signal_duration > 1:
        signal_duration = 1.0    # Max. 1 second
    start = 0.0
    stop = signal_duration     
    step = 1.0 / sample_rate
    t = np.arange(start, stop, step, dtype='float')
    output_buffer = output_buffer[0 : len(t)]
    return (t, output_buffer, len(output_buffer))

def plot_noise(t, s, start_freq, end_freq):
    plt.plot(t, s)
    plt.xlabel('time (s)')
    plt.ylabel('amplitude')
    plt.title('Noise start freq ' + str(start_freq) + ' Hz, stop freq. ' +
              str(end_freq) + ' Hz' )
    plt.grid(True)
    #plt.savefig("test.png")
    plt.show()

def test_noise_signal():
    sample_rate = 40000
    start_freq  = 500.0  
    stop_freq   = 750.0
    signal_duration = 0.01
    t, s, len_s = generate_noise(sample_rate, signal_duration, start_freq, stop_freq)
    plot_noise(t, s, start_freq, stop_freq)
    plot_spectrum(t, s, sample_rate)

########
# FFT Spectrum Analysis.

def plot_spectrum(t, s, sample_rate, ampli_buffer_prev=None, phase_buffer_prev=None):
    if len(s) > sample_rate:
        s = s[0 : sample_rate - 1]

    # realFFT = np.abs(np.fft.rfft(s))
    spectrum = np.fft.rfft(s)
    magnitude = np.abs(spectrum)

    
    for i in range(0, len(magnitude)):
        if magnitude[i] == 0.0:
            magnitude[i] = 1e-20
    output_buffer_magnitude = 10*np.log10(magnitude)

    # Calc magnitude_buffer with prev_buffer in mind.
    if (ampli_buffer_prev is not None):
        buffer_prev_max_value = ampli_buffer_prev.max()
        output_buffer_magnitude = buffer_prev_max_value - ampli_buffer_prev + output_buffer_magnitude

    fft_freq = np.abs(np.fft.fftfreq(len(s), 1 / sample_rate))
    fft_freq = fft_freq[0 : int((len(s) / 2) + 1)]

    # # Artificial limit, only for DEBUG.
    # fft_freq = fft_freq[0 : int(len(fft_freq) / 8)]
    # output_buffer_magnitude = output_buffer_magnitude[0 : int(len(output_buffer_magnitude) / 8)]

    plt.plot(fft_freq, output_buffer_magnitude)
    
    plt.xlabel('freq (Hz)')
    plt.ylabel('amplitude')
    plt.title('FFT Spectrum')
    plt.grid(True)
    #plt.savefig("test.png")
    plt.show()

    output_buffer_phase = np.angle(spectrum)
    
    # Calc phase_buffer with prev_buffer in mind.
    if (phase_buffer_prev is not None):
        buffer_prev_max_value = phase_buffer_prev.max()
        output_buffer_phase = buffer_prev_max_value - phase_buffer_prev + output_buffer_phase
    
    plt.plot(fft_freq, output_buffer_phase)

    plt.ylabel('phase(Radians)')
    plt.title('FFT Phase')
    plt.show()

    ampli_buffer_prev = magnitude
    phase_buffer_prev = output_buffer_phase
    return (ampli_buffer_prev, phase_buffer_prev)

#####################
#####################

# 0. Read configuration file with filter program parameters.

def read_CSV_file_program_parameters():
    filepath = ".\\input_freq_response\\"
    filename = "config.csv"
    
    data = None
    
    try:
        with open(filepath + filename, newline='') as f:
            reader = csv.reader(f)
            data = list(reader)
    except OSError as err:
        print("OS error: {0}".format(err))
    except:
        print("Unexpected error while reading file config.csv:", sys.exc_info()[0])
        exit(0)

    if len(data) < 2:
        print("There are no values in the config.csv file!")
        exit(0)

    data = data[1]

    filepath = data[0]
    filename = data[1]
    num_interpol_possible_values = ["2", "3", "4", "5", "6", "7", "8", "9", "10", "20", "50", "100"]
    if data[2] not in (num_interpol_possible_values):
        print("ERROR in config.csv file, factor_num_of_interpolated_points as to be one of the following values.", num_interpol_possible_values)
        exit(0)
    factor_num_of_interpolated_points = int(data[2])
    flag_possible_values = ["True", "False"]
    if data[3] not in (flag_possible_values):
        print("ERROR in config.csv file, flag_interpolated as to be True or False.")
        exit(0)
    flag_interpolated = eval(data[3])
    N_val = -1
    try:
       N_val = int(data[4])
       if N_val < 0:
           raise ValueError  
    except ValueError:
        print("ERROR in config.csv file, N_val as to be an integer number zero or greater but can be simultaneous zero with D_val.")
        exit(0)
    D_val = -1
    try:
       D_val = int(data[5])
       if D_val < 0 or (D_val == 0 and N_val == 0):
           raise ValueError  
    except ValueError:
        print("ERROR in config.csv file, D_val as to be an integer number zero or greater but can be simultaneous zero with N_val.")
        exit(0)
    sample_rate = -1
    try:
       sample_rate = int(data[6])
       if sample_rate < 10:   # Note: This number could be 2.
           raise ValueError  
    except ValueError:
        print("ERROR in config.csv file, sample_rate as to be an integer number greater or equal to 10.")
        exit(0)
    
    return [ filepath,filename,factor_num_of_interpolated_points,flag_interpolated,N_val,D_val,sample_rate ]
    
#
# 1. Start with a file containing the frequency response that we would like that our filter to have. This CSV file should have many data points. The CSV file has the following structure: <br>
# frequency, amplitude (in linear scale), phase (in radians) <br>
# An interpolation function could be optionally added to generate more intermediate points. <br>
# In here, we will have to choose the M, the number of data points (tuple) that we will be using. The M will dictate the number of rows of the matrix X.
#

def read_CSV_file_frequency_response(filepath, filename):
    data = None

    try:
        with open(filepath + filename, newline='') as f:
            reader = csv.reader(f)
            data = list(reader)
    except OSError as err:
        print("OS error: {0}".format(err))
    except:
        print("Unexpected error while reading file of frequency response parameterized in file config.csv:", sys.exc_info()[0])
        exit(0)

    data = [ [float(point_str[0]), float(point_str[1]), float(point_str[2])] for point_str in data[1:] ]
    return data

def calculate_volt_factor(dBV):
    # 1 V = 0 dBv.
    #
    # The formula for Volts to dBv conversion is:
    #
    # dBV = 20 * log10( Volts )
    #
    # reverse formula for converting dBv to Volts is:
    #
    # Volts = 10 ^ ( dBV/20 )

    volt_scale_factor = math.pow(10.0, dBV / 20.0)
    return volt_scale_factor

def getInterpolated_data_for_freq(data, freq, num_points):
    x_freq_points = np.array([point[0] for point in data])
    y_attenu_points = np.array([point[1] for point in data])
    y_phase_points = np.array([point[2] for point in data])

    # Perform the interpolation.
    interp_attenu = np.interp(freq, x_freq_points, y_attenu_points)
    # TODO: Discomment after DEBUG.
    # NOTE: The attenuation table has to be changed to dBV. 
    #
    # interp_volts_factor = calculate_volt_factor(interp_attenu)
    interp_phase = np.interp(freq, x_freq_points, y_phase_points)
    
    # return (freq, interp_attenu, interp_volts_factor, interp_phase)
    data_interpol = [ [ point[0], point[1], point[2] ] for point in zip(freq, interp_attenu, interp_phase) ]
    return data_interpol


if __name__ == "__main__":
    # test_sine_signal()
    # test_2_sines_signal()
    # test_N_sines_signal()
    # test_chirp_signal()
    # test_dirac_impulse_signal()
    # test_step_signal()
    # test_square_signal()
    # test_square_wave_signal()
    # test_noise_signal()
    

# 0. Read configuration file with filter program parameters.

    parameters = read_CSV_file_program_parameters()
    filepath, filename, factor_num_of_interpolated_points, flag_interpolated, N_val, D_val,sample_rate = parameters
    print("\nParameters from config.csv file: \nfilepath_freq_response: ", filepath,
          "\nfilename_freq_response: ", filename,
          "\nfactor_num_of_interpolated_points: ", factor_num_of_interpolated_points,
          "\nflag_interpolated: ", flag_interpolated,
          "\nN_val: ", N_val,
          "\nD_val: ", D_val,
          "\nsample_rate: ", sample_rate)
    print("\n\n")


    #######
    # 1. Start with a file containing the frequency response that we would like
    # that our filter to have. This CSV file should have many data points. The
    # CSV file has the following structure:
    # frequency, amplitude (in linear scale), phase (in radians)
    # An interpolation function could be optionally added to generate more
    # intermediate points.
    # In here, we will have to choose the M, the number of data points (tuple)
    # that we will be using. The M will dictate the number of rows of the
    # matrix X.

    #######
    # Read transfer function from CSV file.

    # filepath = ".\\input_freq_response\\"
    # filename = "freq_response.csv"
    # filename = "freq_response_book.csv"
    freq_response = read_CSV_file_frequency_response(filepath, filename)

    print(freq_response)

    #######
    # Make interpolation on the data.

    # factor_num_of_interpolated_points = 10    # 10 times the CSV number of points.

    # flag_interpolated = False   # Note: For the book example is False.

    M_val = None

    if (flag_interpolated == True):
        start = 0.0
        stop  = float(freq_response[-1][0])
        num_points = len(freq_response) * 10
        freq = np.linspace(start, stop, num_points)
        freq_resp_interpol = getInterpolated_data_for_freq(freq_response, freq, num_points)

        for point in freq_resp_interpol:
            print(point)

        # We set the current freq_response list to the interpolated version.         
        freq_response = freq_resp_interpol

    M_val = len(freq_response)


    #######
    # 2. Choose the N=? and D=? values. The N is the degree of the polynomial in
    #    the numerator and the D is the degree of the polynomial in the
    #    denominator of equation (1). Only one N or D can be zero at the same time,
    #    but they can be a zero or a positive number of any value. The greater the N
    #    and D number, the larger the filter will be that we are constructing and the
    #    more demanding in terms of computational cost it will be to run it.

    # N_val = 2
    # D_val = 2

    # Sanity check!
    if (N_val == 0 and D_val == 0):
        print("ERROR: Both N and D cannot be zero at the same time!")
        exit(0)


    #######
    # 3. Choose the sample_rate of the digital signal to choose the T_s, that
    #    is the period of the sample rate that is calculated with equation (9).

    # sample_rate =  40000   # Samples / sec
    
    # sample_rate_book = 1000  # Samples / sec
    # sample_rate = sample_rate_book # Samples / sec

    # Period
    T_s = 1 / float(sample_rate)


    #######
    # 4. Fill the y_m sequence of cosines for the length of D.

    def calc__y_m(D_val, T_s, freq_M, A_M, phase_M):
        # D_val   - Polynomial degree of the denominator of the transfer function. 
        # T_s     - Period, inverse of the sample frequency.
        # freq_M  - Frequency of the specific data point.
        # A_M     - Amplitude of the signal at a data point fo a specific frequency.
        # phase_M - Phase of the data point at at a specific frequency. 

        # in here -1 is the D_val polynomial degree term number.
        #    -y_1( -1 ) = A_1 * cos((-1)*W_1*T_s + phi_1) = -A_1 * cos(-W_1*T_s + phi_1)

        y_m_list = []
        for term_num in range(1, D_val + 1):
            W_M = 2*np.pi*freq_M 
            y_m_list.append( -A_M * np.cos( (-term_num) * W_M * T_s + phase_M) )
        return y_m_list


    #######
    # 5. Fill the u_m sequence of cosines for the length of N+1.
    #    (Note: In here is plus one terms).

    def calc__u_m(N_val, T_s, freq_M):
        # N_val   - Polynomial degree of the numerator of the transfer function. 
        # T_s     - Period, inverse of the sample frequency.
        # freq_M  - Frequency of the specific data point.
        
        # in here 0 is the N_val polynomial degree term number.
        #    u_1( -0 ) = cos((-0)*W_1*T_s) = cos(-0*W_1*T_s)

        u_m_list = []
        for term_num in range(0, N_val + 1):
            W_M = 2*np.pi*freq_M 
            u_m_list.append( np.cos( (-term_num) * W_M * T_s) )
        return u_m_list


    #######
    # 6. Join the columns of y_m concatenated with u_m to make the matrix X of
    #    matrix equation X (2).

    # Create the X_matrix filled with zeros. 
    rows = len(freq_response)
    cols = D_val + N_val + 1
    X_matrix = np.zeros((rows, cols))

    for row, data_point in zip(range(0, len(freq_response)), freq_response):
        freq, attenu, phase_shift = data_point
        freq_M  = freq
        A_M     = attenu
        phase_M = phase_shift
        y_m_list = calc__y_m(D_val, T_s, freq_M, A_M, phase_M)

        for y_m, value in zip(range(0, D_val), y_m_list):
            X_matrix[row][y_m] = value
        
        u_m_list = calc__u_m(N_val, T_s, freq_M)
        for u_m, value in zip(range(0, N_val + 1), u_m_list):
            X_matrix[row][D_val + u_m] = value

    print("\n\n X_matrix \n\n", X_matrix)


    #######
    # 7. Create the column vector Y of equation (3).

    # Create the Y_matrix filled with zeros. 
    rows = len(freq_response)
    cols = 1
    Y_matrix = np.zeros((rows, cols))

    #      ---              --- 
    #      | A_1 * cos(phi_1) |
    #  Y = | A_2 * cos(phi_2) |
    #      |      ......      |
    #      | A_m * cos(phi_m) |
    #      ---              ---

    for row, data_point in zip(range(0, len(freq_response)), freq_response):
        freq, attenu, phase_shift = data_point
        A_M   = attenu
        phi_M = phase_shift
        Y_matrix[row][0] = A_M * np.cos(phi_M)

    print("\n\n Y_matrix \n\n", Y_matrix)


    #######
    # 8. Calculate the Pseudo-Inverse the vector THETA. It will have the filter
    #    coefficients. Depending on the number chosen for N and D, a different
    #    type of filter topology will be made.


    THETA = pinv(X_matrix).dot( Y_matrix)    
    # THETA = inv(np.transpose(X_matrix) * X_matrix) * np.transpose(X_matrix) * Y_matrix    

    print("\n\n THETA matrix: \n\n", THETA)

    #######
    # 9. Construct the filter.

    # Get the a_D and b_N filter coefficients. 
    a_0__start = [1.0]
    a_D_filter_coef = a_0__start + [elem[0] for elem in THETA[0:D_val].tolist()] 
    b_N_filter_coef = [elem[0] for elem in THETA[D_val:].tolist()]

    print("\n\na_D__filter_coefficients: ")
    for i, a_val in enumerate(a_D_filter_coef):
        print("  a_{0}: {1:.4f}".format(i, a_val))

    print("\n\nb_N__filter_coefficients: ")
    for i, b_val in enumerate(b_N_filter_coef):
        print("  b_{0}: {1:.4f}".format(i, b_val))


    def our_filter(a_D_coef, b_N_coef, in_sequence):
        out_sequence = []
        a0, a1, a2 = a_D_coef
        b0, b1, b2 = b_N_coef
        p0 = 0.0
        p1 = 0.0
        p2 = 0.0
        for elem in in_sequence[-1::-1]:
            p0  = elem + (-a1*p1) + (-a2*p2)
            out = b0*p0 + b1*p1 + b2*p2  
            p2 = p1
            p1 = p0
            out_sequence.insert(0, out)
        return out_sequence


    def our_filter_generic(a_D_coef, b_N_coef, in_sequence, generate_C_filter_file = False):
        # Makes a generic filter from the list of coefficients.

        out_sequence = []

        # Not generic.
        # a0, a1, a2 = a_D_coef
        # b0, b1, b2 = b_N_coef
        # p0 = 0.0
        # p1 = 0.0
        # p2 = 0.0

        # Generic.
        a_D = a_D_coef
        b_N = b_N_coef
        p_T = [0.0 for elem in max(a_D_coef, b_N_coef)] 
        
        for elem in in_sequence[-1::-1]:        
            # Not generic.
            # p0  = elem + (-a1*p1) + (-a2*p2)
            # out = b0*p0 + b1*p1 + b2*p2  
            # p2 = p1
            # p1 = p0

            # Generic.
            p_T[0] = elem + sum( [ -a*p for a, p in zip( a_D[1:len(a_D)], p_T[1:len(a_D)] ) ] )
            out = sum( [ b*p for b, p in zip( b_N[0:len(b_N)], p_T[0:len(b_N)] ) ] )
            p_T[1:] = p_T[0:-1]

            out_sequence.insert(0, out)

        # Generate the C output file.
        if generate_C_filter_file == True:
            write_src_code_to_file(a_D_coef, b_N_coef)

        return out_sequence

    def write_src_code_to_file(a_D_coef, b_N_coef):
        # Generate the C output file.
        filepath = ".\\"
        dot_H_file_name = "custom_filter.h"
        src_code = generate_C_source_code_filter_file(a_D_coef, b_N_coef, dot_H_file_name)
        print("\n\n")
        print(src_code)
        with open(filepath + dot_H_file_name, "w") as f:
            f.write(src_code)
        return src_code

    def generate_C_source_code_filter_file(a_D_coef, b_N_coef, dot_H_file_name):
        # Makes a generic C filter file from the list of coefficients.

        a_D = a_D_coef
        b_N = b_N_coef

        # .H guard
        src_code =  "#ifndef %s_H_\n"    % (dot_H_file_name[:-2].upper())
        src_code += "#define %s_H_\n\n"  % (dot_H_file_name[:-2].upper())
        
        src_code += '''//generated_C_.h_source_code_filter_file\n\n\n''' 

        # Constants and variables.
        for i, coef in zip(range(1, len(a_D_coef)), a_D_coef[1:]):   # enumerate(a_D_coef):
            src_code += "float const a_%d = %f;\n" % (i, coef)  
        src_code += "\n"

        for i, coef in enumerate(b_N_coef):
            src_code += "float const b_%d = %f;\n" % (i, coef)  
        src_code += "\n"
        
        for i, _ in enumerate(max(a_D_coef, b_N_coef)):
            src_code += "float p_%d = (float) 0.0;\n" % (i)  
        src_code += "\n"

        # Function custom_c_filter_init() 
        src_code += "void custom_filter_init(){\n"

        for i, _ in enumerate(max(a_D_coef, b_N_coef)):
            src_code += "    float p_%d = 0.0;\n" % (i)  
        # src_code += "\n"

        src_code += "}\n\n"

        src_code += "float custom_c_filter(float inputValue){\n"

        # p_T[0] = elem + sum( [ -a*p for a, p in zip( a_D[1:len(a_D)], p_T[1:len(a_D)] ) ] )
        src_code += "    p_0 = inputValue "
        for a, p in zip( range(1, len(a_D)), range(1, len(a_D)) ):
            src_code += " + (-a_%d*p_%d)" % (a, p)
        src_code += ";\n"

        # out = sum( [ b*p for b, p in zip( b_N[0:len(b_N)], p_T[0:len(b_N)] ) ] )
        src_code += "    float outputValue = "
        for b, p in zip( range(0, len(b_N)), range(0, len(b_N)) ):
            src_code += " + (b_%d*p_%d)" % (b, p)
        src_code += ";\n"
        
        # p_T[1:] = p_T[0:-1]
        for i in range(len(max(a_D_coef, b_N_coef)[:-1]) - 1, -1, -1):
            src_code += "    p_%d = p_%d;\n" % (i+1, i)  
        # src_code += "\n"
            
        src_code += "    return outputValue;\n}\n\n"
        
        # .H guard
        src_code += "#endif\n"
        
        return src_code


    #######
    # 10. Test the filter frequency response with a simples signal impulse. 
        
    in_sequence = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0]
    # out_sequence = our_filter(a_D_filter_coef, b_N_filter_coef, in_sequence)
    generate_C_filter_file = True
    out_sequence = our_filter_generic(a_D_filter_coef, b_N_filter_coef, in_sequence, generate_C_filter_file)
    generate_C_filter_file = False

    print("\n\nImpulse signal response")
    print("in_sequence: ",  in_sequence)
    print("out_sequence: ", out_sequence)
    print("\n")


    #######
    # 11. Test the filter frequency response with a chirp signal between the
    #     frequency limits.

    # Generate Chirp.
    sample_rate    = 1000
    start_freq     = 0.0
    end_freq       = 500.0
    chirp_duration = 1.00
    t, s, len_s = generate_chirp(sample_rate, chirp_duration, start_freq, end_freq)
    
    # plot_chirp(t, s, start_freq, end_freq)
    # plot_spectrum(t, s, sample_rate)    

    in_sequence = s.tolist()
    # out_sequence = our_filter(a_D_filter_coef, b_N_filter_coef, in_sequence)
    out_sequence = our_filter_generic(a_D_filter_coef, b_N_filter_coef, in_sequence, generate_C_filter_file)

    s_new = np.array(out_sequence)

    # plot_chirp(t, s_new, start_freq, end_freq)
    # plot_spectrum(t, s_new, sample_rate)    


    #######
    # 12. Test the filter frequency response with a N_sines_signal() signal between the
    #     frequency limits.

    sample_rate     = 1000
    freq_list       = [float(freq) for freq in range(0, int(sample_rate / 2.0), 10)]
    phase_list      = [0.0  for freq in range(0, int(sample_rate / 2.0), 10)]
    amplitude_list  = [1.0  for freq in range(0, int(sample_rate / 2.0), 10)]
    signal_duration = 10
    t, s, len_s = generate_N_sines(sample_rate, signal_duration, freq_list, phase_list, amplitude_list)
    
    plot_N_sines(t, s, freq_list, phase_list)
    ampli_buffer_prev, phase_buffer_prev = plot_spectrum(t, s, sample_rate)

    in_sequence = s.tolist()
    # out_sequence = our_filter(a_D_filter_coef, b_N_filter_coef, in_sequence)
    out_sequence = our_filter_generic(a_D_filter_coef, b_N_filter_coef, in_sequence, generate_C_filter_file)

    s_new = np.array(out_sequence)

    plot_chirp(t, s_new, start_freq, end_freq)
    plot_spectrum(t, s_new, sample_rate, ampli_buffer_prev, phase_buffer_prev)


    #######
    # 13. Test the filter frequency response with a noise_signal() signal between the
    #     frequency limits.

    sample_rate     = 1000
    start_freq      = 0.0  
    stop_freq       = 500.0
    signal_duration = 1
    t, s, len_s = generate_noise(sample_rate, signal_duration, start_freq, stop_freq)
    
    # plot_noise(t, s, start_freq, stop_freq)
    # ampli_buffer_prev, phase_buffer_prev = plot_spectrum(t, s, sample_rate)

 
    in_sequence = s.tolist()
    # out_sequence = our_filter(a_D_filter_coef, b_N_filter_coef, in_sequence)
    out_sequence = our_filter_generic(a_D_filter_coef, b_N_filter_coef, in_sequence, generate_C_filter_file)

    s_new = np.array(out_sequence)

    # plot_chirp(t, s_new, start_freq, end_freq)
    # plot_spectrum(t, s_new, sample_rate, ampli_buffer_prev, phase_buffer_prev)    


