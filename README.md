# Radar-Target-Generation-Detection
Fourth Project of SFND

## Project Layout

<img src="/images/layout.png" width="700" />

1. Configure the FMCW waveform based on the system requirements.
2. Define the range and velocity of target and simulate its displacement.
3. For the same simulation loop process the transmit and receive signal to determine the beat signal
4. Perform Range FFT on the received signal to determine the Range
5. Towards the end, perform the CFAR processing on the output of 2nd FFT to display the target.


## 1. Radar Specifications

- Frequency of operation = 77GHz
- Max Range = 200m
- Range Resolution = 1 m
- Max Velocity = 100 m/s

```
max_range = 200;
range_res = 1;
max_vel = 100;
vel_res = 3;
C = 3e8; 
```

## 2. User Defined Range and Velocity of target
- define the target's initial position and velocity. 
- Note : Velocity remains contant
```
v_0 = 30; % 30 m/s
target_pos = 100; % 100 m 
```

## 3. FMCW Waveform Generation

- Design the FMCW waveform by giving the specs of each of its parameters.
- Calculate the Bandwidth (B), Chirp Time (Tchirp) and Slope (slope) of the FMCW chirp using the requirements above.
```
b_sweep = C / (2*range_res);
Ts = (2*max_range) / C;
Ts = 5.5*Ts; % For an FMCW radar system, the sweep time should be at least 5 to 6 
             % times the round trip time.
slope = b_sweep / Ts;

```

- The number of chirps in one sequence. Its ideal to have 2^ value for the ease of running the FFT for Doppler Estimation. 
```
Nd=128;                   % #of doppler cells OR #of sent periods % number of chirps
```

- The number of samples on each chirp. 
```
Nr=1024;                  %for length of time OR # of range cells
```

- Timestamp for running the displacement scenario for every sample on each chirp
```
t=linspace(0,Nd*Tchirp,Nr*Nd); %total time for samples
```

## 4. Signal generation and Moving Target simulation

- Running the radar scenario over the time. 
```
for i=1:length(t)         
    
    
    % *%TODO* :
    %For each time stamp update the Range of the Target for constant velocity. 
    r_t(i) = target_pos + v_0 * t(i); 
    td(i) = 2*r_t(i) / C;
    
    % *%TODO* :
    %For each time sample we need update the transmitted and
    %received signal. 
    Tx(i) = cos(2*pi*(fc*t(i) + 0.5*slope*t(i)^2));
    Rx(i) = cos(2*pi*(fc*(t(i) - td(i)) + 0.5*slope*(t(i)-td(i))^2 ));
     
    % *%TODO* :
    %Now by mixing the Transmit and Receive generate the beat signal
    %This is done by element wise matrix multiplication of Transmit and
    %Receiver Signal
    
    Mix(i) = Tx(i) .* Rx(i);
    
end
```

## 5. RANGE MEASUREMENT

- reshape the vector into Nr*Nd array. Nr and Nd here would also define the size of Range and Doppler FFT respectively.
```
mix_sig_reshaped = reshape(Mix, Nr, Nd);
```

- run the FFT on the beat signal along the range bins dimension (Nr) and normalize.
```
sig_fft_1 = fft(Mix_reshape, Nr);
sig_fft_1 = sig_fft_1 ./ Nr;
```

- Take the absolute value of FFT output
```
sig_fft_1 = abs(sig_fft_1);
```

- Output of FFT is double sided signal, but we are interested in only one side of the spectrum. Hence we throw out half of the samples.
```
sig_fft_1 = sig_fft_1(1:Nr/2);
```

## 6. RANGE DOPPLER RESPONSE

- The 2D FFT implementation is already provided here. 
- This will run a 2DFFT on the mixed signal (beat signal) output and generate a range doppler map.
- You will implement CFAR on the generated RDM Range Doppler Map Generation.
- The output of the 2D FFT is an image that has reponse in the range and doppler FFT bins. 
- So, it is important to convert the axis from bin sizes to range and doppler based on their Max values.
```
Mix = reshape(Mix,[Nr,Nd]);
```
* 2D FFT using the FFT size for both dimensions.
```
sig_fft2 = fft2(Mix,Nr,Nd);
```
* Taking just one side of signal from Range dimension.
```
sig_fft2 = sig_fft2(1:Nr/2,1:Nd);
sig_fft2 = fftshift (sig_fft2);
RDM = abs(sig_fft2);
RDM = 10*log10(RDM) ;
```
* Use the surf function to plot the output of 2DFFT and to show axis in both dimensions
```
doppler_axis = linspace(-100,100,Nd);
range_axis = linspace(-200,200,Nr/2)*((Nr/2)/400);
figure ('Name','Range and Speed From FFT2')
surf(doppler_axis,range_axis,RDM);
```

## 7. CFAR implementation

- Slide Window through the complete Range Doppler Map
- Select the number of Training Cells in both the dimensions.
```
Tr = 14;
Td = 8;
```


-Select the number of Guard Cells in both dimensions around the Cell under test (CUT) for accurate estimation
```
Gr = 3;
Gd = 2;
```


- offset the threshold by SNR value in dB
```
thre_offset = 1.1;
```

- a loop such that it slides the CUT across range doppler map by giving margins at the edges for Training and Guard Cells.
- For every iteration sum the signal level within all the training cells. To sum convert the value from logarithmic to linear using db2pow function. 
- Average the summed values for all of the training cells used. After averaging convert it back to logarithimic using pow2db.
- Further add the offset to it to determine the threshold. 
- Next, compare the signal under CUT with this threshold. If the CUT level > threshold assign it a value of 1, else equate it to 0.
- Use RDM[x,y] as the matrix from the output of 2D FFT for implementing
 
```
RDM = RDM/max(max(RDM));

for i = Tr+Gr+1:(Nr/2)-(Gr+Tr)
    for j = Td+Gd+1:Nd-(Gd+Td)
       
        % vector to store noise_level for each iteration on training cells
        noise_level = zeros(1,1);
        
        for p = i-(Tr+Gr) : i+(Tr+Gr)
            for q = j-(Td+Gd) : j+(Td+Gd)
                if (abs(i-p) > Gr || abs(j-q) > Gd)
                    noise_level = noise_level + db2pow(RDM(p,q));
                end
            end
        end
        
        threshold = pow2db(noise_level/(2*(Td+Gd+1)*2*(Tr+Gr+1)-(Gr*Gd)-1));
        threshold = threshold + thre_offset;
        CUT = RDM(i,j);
        
        if (CUT < threshold)
            RDM(i,j) = 0;
        else
            RDM(i,j) = 1;
        end
        
    end
end
```

- The process above will generate a thresholded block, which is smaller  than the Range Doppler Map as the CUT cannot be located at the edges ofmatrix. Hence,few cells will not be thresholded. To keep the map size same set those values to 0. 
```
[lr, lc] = size(RDM);
RDM(union(1:(Tr+Gr), lr-(Tr+Gr-1):lr), :) = 0;  
RDM(:, union(1:(Td+Gd), lc-(Td+Gd-1):lc)) = 0;  
```


- display the CFAR output using the Surf function like we did for Range Doppler Response output.
```
figure,surf(doppler_axis,range_axis,RDM);
colorbar;
```

