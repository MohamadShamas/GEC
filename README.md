# GEC
Gamma Event Coupling code + dataset

In orfder to run GEC, user must first read the EDF files using the povided fileio functions (freely availble with Fieldtrip toolbox). The usage syntax is as follows:

header = ft_read_header(filename);
Sig = ft_read_data(filename,'heasder', header);

Then, sampling frequency can be extracted from the header structure:

Fs = header.Fs

Finally, GEC matrices can be computed using the following syntax:

sWindow = 30 % window size in seconds
sStep = 30 % step size in seconds (here 0 overlap)
Fmin = 30 % lower frequncy band limit, in case of low gamma it is 30HZ
Fmax = 55 % upper frequncy band limit, in case of low gamma it is 55Hz

gec = GEC(Sig,sWind,sStep,Fs,Fmin,Fmax)
