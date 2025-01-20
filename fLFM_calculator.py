import numpy as np

print('## The purpose of this calculator is to help design a fourier light field microscope.')
print('## The calculations are based on Galdón et al., 2022. Fourier lightfield microscopy: a practical design guide. DOI: 10.1364/ao.453723')
print('* Start with the imaging specs you want to achieve (hard-coded in the script)')
print('* Decide on a camera')
print('* Find a microlens arrray/lenslet array that fits')
print('* Choose objective lens and aperture stop')
print('* Calculate relay and field stop\n')

#functions
def pitch_suggestion(sensor_size,N):
    P = sensor_size/N
    return P

def NAef_of_MLA(P,fMLA):
    NAef = np.sin(np.arctan(P/fMLA/2))
    return NAef

def necessary_total_magnification_for_FOV(FOV_target,P):
    Mt = P/FOV_target
    return Mt

def necessary_total_magnification_for_resolution(Lambda, resolution_target,NAef,delta):
    Mt = 2*delta/(resolution_target-Lambda/NAef/2)
    return Mt

def necessary_total_magnification_for_DOF(Lambda, DOF_target,NAef,delta):
    Mt = delta*NAef/(DOF_target*NAef**2 - Lambda)
    return Mt

def depth_of_field(Lambda_emission,NAef,delta,Mt):
    DOF = 2*Lambda_emission/(NAef**2) + delta/(Mt*NAef)
    return DOF

def NAef_given_target_resolution(Lambda_emission,resolution_target,Mt,delta):
    NAef = Lambda_emission/(2*(resolution_target-2*delta/Mt))
    return NAef

def Nef_given_target_DOF(Lambda_emission,DOFef_target,Mt,delta):
    NAef = (np.sqrt(delta**2 + 4*DOFef_target*Lambda_emission*Mt**2) + delta)/(2*DOF_target*Mt)
    return NAef

def resolution(Lambda_emission,NAef,delta,Mt):
    res = Lambda_emission/(2*NAef) + 2*delta/Mt
    return res

def relay_magnification(Mt,fMLA,fMO):
    Mr = fMLA/(fMO*Mt)
    return Mr
    
def second_relay_f(beam_path,Mr):
    f2 = Mr*beam_path/(2*Mr+2)
    return f2

def first_relay_f(f2,Mr):
    f1 = f2/Mr
    return f1

def field_stop(P,f2,fMLA):
    FS = P*f2/fMLA
    return FS

def aperture_stop(NAobj,fMO):
    ASobj = 2*np.tan(np.arctan(NAobj))*fMO
    return ASobj

def field_of_view(P,Mt):
    FOV = P/Mt
    return FOV

def paralax_angle(f1,f2,fMO,P):
    sigma = P*f1/f2/fMO
    return np.rad2deg(sigma)

def delta_z(fMLA, P, delta, Mt):
    delta_z = (fMLA/P)*delta/Mt**2
    return delta_z

# LETS GO!!
print('Targets and constraints')
Lambda_emission = 0.00051 # in nm, 510nm
size_of_neuron = 0.015 # in mm, 15um
sampling_factor = 2.5 # choose in pixels per cell, just a little higher than the Nyquist rate
resolution_target_suggestion = size_of_neuron/sampling_factor# in mm, `6um is alternative
FPV_target_suggestion = 3 # in mm, 3mm
DOF_target_suggestion = 0.2 # in mm, 200um
N=float(input('Number of elemental images (5, 3 or 7): ').strip() or '5') # INPUT

print('\n')

print('decide on camera specs')
Mpix = float(input('Camera megapixels (eg 67): ').strip() or '67') # INPUT
delta = float(input('Pixel size in um (eg 2.5): ').strip() or '2.5') # INPUT
delta = delta/1000 # in mm
sensor_size =  np.sqrt(Mpix*10**6)*delta # in mm
print('\nCamera: ' + str(Mpix) + ' Mpix, ' + str(delta) + ' um pixel size, ' + str(round(sensor_size,2)) + ' mm sensor size')
print('\n')

# choose MLA
# calculations recommended MLA pitch
P = pitch_suggestion(sensor_size,N)
print('maximum MLA pitch = ' + str(round(P,3)) + ' mm')
print('find MLA with pitch as close to max as possible and enter pitch and fMLA')
P = float(input('Choose P in mm (<' + str(round(P,4)) +'): ').strip() or str(round(P,4))) # INPUTS
fMLA = float(input('Choose fMLA in mm (~10xP): ').strip() or str(round(10*P,2))) # INPUTS
NAef = NAef_of_MLA(P,fMLA)
print('MLA pitch = ' + str(round(P,2)) + ' mm, focal length = ' + str(fMLA) + ' mm, NA = ' + str(round(NAef,2)))
print('min objective NA to match MLA NA: ' + str(round(NAef*N,2)))

# choose objective 
fMO = float(input('Choose fMO in mm (for Nikon 10X/0.3NA = 20): ').strip() or '20') #INPUTS
ASobj_min = aperture_stop(NAef*N,fMO)


choose_target_parameter = input('Choose target resolution (r), target field of view (f), \n,' + 
                                'target depth of field (d), or calculate specs based on given relay (c): ').strip() or 'r'
if choose_target_parameter == 'r':
    resolution_target = float(input('target resolution in um (eg. ' + str(round(resolution_target_suggestion*1000,1)) + '):').strip() or 
                          str(round(resolution_target_suggestion*1000,1)))/1000 # in mm 
    Mt = necessary_total_magnification_for_resolution(Lambda_emission,resolution_target,NAef,delta)
    FOV = field_of_view(P,Mt)
    DOF = depth_of_field(Lambda_emission,NAef,delta,Mt)
    print('necessary total magnification to get resolution = x' + str(round(Mt,2)))
    print('field of view = ' + str(round(FOV,2)) + ' mm')
    print('depth of field = ' + str(round(DOF,2)) + ' mm')

elif choose_target_parameter == 'f':
    FOV_target = float(input('target field of view in mm (eg. 3): ').strip() or '3')
    Mt = necessary_total_magnification_for_FOV(FOV_target,P)
    res = resolution(Lambda_emission,NAef,delta,Mt)
    DOF = depth_of_field(Lambda_emission,NAef,delta,Mt)
    print('necessary total magnification to get FOV = x' + str(round(Mt,2)))
    print('resolution = ' + str(round(res*1000,2)) + ' um')
    print('depth of field = ' + str(round(DOF,2)) + ' mm')

elif choose_target_parameter == 'd':
    DOF_target = float(input('target depth of field in mm (eg. 0.2): ').strip() or '0.2')
    Mt = necessary_total_magnification_for_DOF(Lambda_emission, DOF_target,NAef,delta)
    res = resolution(Lambda_emission,NAef,delta,Mt)
    DOF = depth_of_field(Lambda_emission,NAef,delta,Mt)
    print('necessary total magnification to get FOV = x' + str(round(Mt,2)))
    print('resolution = ' + str(round(res*1000,2)) + ' um')
    print('depth of field = ' + str(round(DOF,2)) + ' mm')

# if  choose_target_parameter not c
if choose_target_parameter != 'c':
    # determine relay
    beam_path = float(input('\nChoose beam path length in mm (~600): ').strip() or '600')
    Mr = relay_magnification(Mt,fMLA,fMO)
    f2 = second_relay_f(beam_path,Mr)
    f1 = first_relay_f(f2,Mr)
    FS = field_stop(P,f2,fMLA)
elif choose_target_parameter == 'c':
    f1 = float(input('f1 in mm (eg 200): ').strip() or '200')
    f2 = float(input('f2 in mm (eg 100): ').strip() or '100')
    Mr = f2/f1
    FS = field_stop(P,f2,fMLA)
    Mt = fMLA/fMO/Mr
print('\n')

# recalculate microscope specs
res = resolution(Lambda_emission,NAef,delta,Mt)
DOF = depth_of_field(Lambda_emission,NAef,delta,Mt)
FOV = field_of_view(P,Mt)
sigma = paralax_angle(f1,f2,fMO,P)
delta_z = delta_z(fMLA, P, delta, Mt)

# relay lens diameters
min_diam_f1 = FOV*f1/fMO # mm (not sure this is correct)
min_diam_f2 = min_diam_f1*f2/f1 # mm (not sure this is correct)

# summary
print('\n----- Final Specs ----')
print('\nComponents')
print('Camera: ' + str(Mpix) + ' Mpix, pixel size ' + str(delta) + ' mm, sensor size ' + str(round(sensor_size,2)) + ' mm')
print('MLA: pitch = ' + str(round(P,3)) + ' mm, focal length = ' + str(round(fMLA,2)) + ' mm', 'NA = ' + str(round(NAef,2)))
print('Objective: focal length ' + str(round(fMO,2)) + ' mm, objective aperture stop = ' + str(round(ASobj_min,2))
       + ' mm, min NA = ' + str(round(NAef*N,2)))
print('Relay: f1 = ' + str(round(f1,2)) + ' mm (min diam ' + str(round(min_diam_f1)) + ' mm),' +
    'f2 = ' + str(round(f2,2)) + ' mm(min diam ' + str(round(min_diam_f2)) + ' mm),' +
    'field stop = ' + str(round(FS,2)) + ' mm')

print('\nSpecs')
print('number of elemental images = ' + str(N))
print('Resolution = ' + str(round(res*1000,2)) + ' um')
print('um per pixel = ' + str(round(delta/Mt*1000,2)) + ' um')
print('depth of field = ' + str(round(DOF,2)) + ' mm')
print('field of view = ' + str(round(FOV,2)) + ' mm')
print('total magnification = x' + str(round(Mt,2)), 'relay magnification = ' + str(round(Mr,2)))
print('effective NA = ' + str(round(NAef,2)))
print('parallax angle = ' + str(round(sigma,2)) + ' degrees')
print('axial distance related to 1-pixel shift = ' + str(round(delta_z*1000,2)) + ' um')
print('\n')
print('------------------------')
print('\n')
