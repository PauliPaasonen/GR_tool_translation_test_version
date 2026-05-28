

function events = growth_events(data)
% GROWTH_EVENTS Detect new particle formation events
%
% INPUT:
%   data.time  - datetime vector
%   data.size  - particle size bins (nm)
%   data.N     - concentration matrix (size x time)
%
% OUTPUT:
%   events - struct array with detected events

time  = data.time;
sizev = data.size;
N     = data.N;

nt = length(time);
ns = length(sizev);

% -----------------------------
% PARAMETERS (from Python logic assumption)
% -----------------------------
min_size = 3;       % nm
max_size = 100;     % nm region of interest
threshold = 1e3;    % concentration threshold

min_event_duration = 3; % number of time steps (≈1.5 h for 30 min data)

% -----------------------------
% FIND SIZE RANGE INDICES
% -----------------------------
idx = find(sizev >= min_size & sizev <= max_size);

if isempty(idx)
    error('No size bins in selected range');
end

% Extract relevant size range
Nsel = N(idx, :);
sizes_sel = sizev(idx);

% -----------------------------
% EVENT DETECTION
% -----------------------------
is_event = false(1, nt);

for t = 1:nt
    profile = Nsel(:, t);

    % condition: concentration above threshold
    if max(profile) > threshold
        is_event(t) = true;
    end
end

% -----------------------------
% GROUP INTO CONTINUOUS EVENTS
% -----------------------------
events = struct();
event_id = 0;

t = 1;

while t <= nt
    
    if is_event(t)
        t_start = t;
        
        % extend until event ends
        while t <= nt && is_event(t)
            t = t + 1;
        end
        
        t_end = t - 1;
        
        duration = t_end - t_start + 1;
        
        % only accept sufficiently long events
        if duration >= min_event_duration
            
            event_id = event_id + 1;
            
            events(event_id).time = time(t_start:t_end);
            events(event_id).indices = t_start:t_end;
            events(event_id).data = Nsel(:, t_start:t_end);
            events(event_id).sizes = sizes_sel;
            
        end
        
    else
        t = t + 1;
    end
    
end

% -----------------------------
% POST-PROCESSING (growth check)
% -----------------------------
% Try to detect size growth (increasing peak size over time)

for i = 1:length(events)
    
    ev = events(i);
    nt_ev = length(ev.indices);
    
    peak_sizes = zeros(1, nt_ev);
    
    for t = 1:nt_ev
        profile = ev.data(:, t);
        
        [~, imax] = max(profile);
        peak_sizes(t) = ev.sizes(imax);
    end
    
    events(i).peak_sizes = peak_sizes;
    
    % Check monotonic growth (rough criterion)
    dsize = diff(peak_sizes);
    
    % NOTE:
    % original Python logic may use more advanced filtering here
    events(i).is_growth = sum(dsize > 0) > length(dsize) / 2;
    
end

% -----------------------------
% FILTER ONLY TRUE GROWTH EVENTS
% -----------------------------
keep = [events.is_growth];

events = events(keep);

% -----------------------------
% NOTE FLAGS
% -----------------------------
% NOTE: 
% - Threshold value should be verified against Python code
% - Event detection likely more sophisticated in original (e.g. smoothing)
% - aerosol-functions may influence preprocessing
% - Time resolution assumed 30 min

end
