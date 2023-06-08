% Basic voice morphing script
% Voice Neurocognition Laboratory
% based on functions wrtiten by Hideki Kawahara and Julien Rouger

%% create M-objects from sounds
mObject=ExtractMObject('french_natural_pitch_intensity.wav') 
% opens spectrogram to set anchors 
% right click for temporal, then left for specTRAL; 
% ESC to start again, ENTER to accept
save MObj6_neutral.mat mObject

mObject=ExtractMObject('6_pleasure.wav') 
save MObj_french_natural.mat mObject

%% Interpolate between M-objects
% define mobjects to interpolate (here 2)
mobjs=cell(1,2);
load Mobj6_neutral.mat; mobjs{1,1}=mObject;
load Mobj6_pleasure.mat; mobjs{1,2}=mObject;

% define weight of each mobject in the interporaltion (here 50/50)
rates.F0=[0.5 0.5]; 
rates.spectralamplitude=[0.5 0.5];
rates.aperiodicity=[0.5 0.5];
rates.time=[0.5 0.5];
rates.frequency=[0.5 0.5];

% execute interpolation
mObjectM=voicemultimorph(mobjs,rates);


%% write new mObject and sound
sy=executeSTRAIGHTsynthesisM(mObject); 
sy=.95*sy/max(abs(sy)); % normalise
mObject.waveform=sy;

% save new M object
save Mobj6_Neu50Plea50.mat mObjectM

% play and rite sound
wavplay(sy,48000); % 
audiowrite('French_natural_resynth.wav',sy,44100);
