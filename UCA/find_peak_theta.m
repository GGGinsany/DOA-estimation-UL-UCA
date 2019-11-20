function esti_theta = find_peak_theta( music_spectr, theta_range, guessed_num )

[~,loc_peaks] = findpeaks( abs(music_spectr),'NPeaks',guessed_num,'SortStr','descend');

esti_theta = theta_range( loc_peaks );