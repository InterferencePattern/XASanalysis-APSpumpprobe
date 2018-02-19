function [ detectorindices ] = find_detector_indices( final_detectornames, photoncounting )
%FIND_DETECTOR_INDICES 
%   Takes in the list of names and gives an index for each value.
detectorindices.TimeIndex = find(strcmp(final_detectornames,'  7idd:DG1:ADelayAO, Chan A Delay, TABLE, Sec, 7idd:DG1:ADelayAI, Chan A Delay, Sec'));
detectorindices.PhaseShifter = find(strcmp(final_detectornames,'  7idd:digDelay1:actualDelay.VAL, , ps'));
detectorindices.I0 = find(strcmp(final_detectornames,'  7idd:scaler1.S2, , '));
detectorindices.EnergyIndex = find(strcmp(final_detectornames,'  7ida:BraggEAO.VAL, Kohzu Energy, TABLE, keV, 7ida:BraggERdbkAO, Kohzu Energy Readback, keV'));

if photoncounting
    % Photon counting ("gated") signal:
    detectorindices.TFY_LaserOFF3sec = find(strcmp(final_detectornames,'  amofg1a:scaler1.S9, , ')); % Photon counting ("gated") signal
    detectorindices.TFY_LaserON3sec = find(strcmp(final_detectornames,'  amofg1a:scaler1.S10, , ')); % Photon counting ("gated") signal
    detectorindices.TFY_LaserOFF6sec = find(strcmp(final_detectornames,'  amofg1a:scaler1.S9, , ')); % Photon counting ("gated") signal
    detectorindices.TFY_LaserON6sec = find(strcmp(final_detectornames,'  amofg1a:scaler1.S10, , ')); % Photon counting ("gated") signal
    % Set fallback in case photon counting is missing:
    detectorindices.TFY_LaserOFF_bak = find(strcmp(final_detectornames,'  7idb:userCalc11.VAL, BoxcarOFF (3s), ')); % Full analog signal, not photon counting
    detectorindices.TFY_LaserON_bak = find(strcmp(final_detectornames,'  7idb:userCalc12.VAL, BoxcarON (3s), ')); % Full analog signal, not photon counting
else
    % % Full analog signal, not photon counting
    detectorindices.TFY_LaserOFF3sec = find(strcmp(final_detectornames,'  7idb:userCalc11.VAL, BoxcarOFF (3s), ')); % Full analog signal, not photon counting
    detectorindices.TFY_LaserON3sec = find(strcmp(final_detectornames,'  7idb:userCalc12.VAL, BoxcarON (3s), ')); % Full analog signal, not photon counting
    detectorindices.TFY_LaserOFF6sec = find(strcmp(final_detectornames,'  7idb:userCalc13.VAL, boxcarOFF (6s), ')); % Full analog signal, not photon counting
    detectorindices.TFY_LaserON6sec = find(strcmp(final_detectornames,'  7idb:userCalc14.VAL, boxcarON (6s), ')); % Full analog signal, not photon counting
end

detectorindices.HERFD_LaserOFF3sec = find(strcmp(final_detectornames,'  amofg1a:scaler1.S2, , '));
detectorindices.HERFD_LaserON3sec = find(strcmp(final_detectornames,'  amofg1a:scaler1.S3, , '));
detectorindices.HERFD_LaserOFF6sec = find(strcmp(final_detectornames,'  amofg1a:scaler1.S5, , '));
detectorindices.HERFD_LaserON6sec = find(strcmp(final_detectornames,'  amofg1a:scaler1.S6, , '));

% Choose 3 or 6 second integration
detectorindices.TFY_LaserON = detectorindices.TFY_LaserON3sec; if isempty(detectorindices.TFY_LaserON); error('TFY Laser ON column not found!'); end;
detectorindices.TFY_LaserOFF = detectorindices.TFY_LaserOFF3sec; if isempty(detectorindices.TFY_LaserON); error('TFY Laser OFF column not found!'); end;
detectorindices.HERFD_LaserON = detectorindices.HERFD_LaserON3sec; if isempty(detectorindices.TFY_LaserON); error('HERFD Laser ON column not found!'); end;
detectorindices.HERFD_LaserOFF = detectorindices.HERFD_LaserOFF3sec; if isempty(detectorindices.TFY_LaserON); error('HERFD Laser OFF column not found!'); end;

end

