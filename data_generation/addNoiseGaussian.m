function v_noisy = addNoiseGaussian(v_clean,focal_length,pixel_noise)
    % ADDNOISEGAUSSIAN This function takes the observation
    % v_clean (in pixels) and add noise
    % sample from a gaussian distribution
    % with sigma pixel_noise


    % asyntropic noise for the observation
    noise_sample = normrnd(0,pixel_noise,[2, 1]);

    v_noisy = v_clean ./ v_clean(3);
    v_noisy(1:2) = v_noisy(1:2) + noise_sample;

end
