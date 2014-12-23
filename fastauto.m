function [P, W] = fastauto(s, winsize, jump, corrthresh)

% [P, W] = fastauto(s, winsize, jump, corrthresh)
%
% Implements a pitch detector using a fast-autocorrelation algorithm.
%
% Returns P, a vector of frequencies corresponding to each point in time
% and W, a vector of frequencies corresponding to each window
%
% Takes s, the monaural signal; winsize, the length of the window; jump, the
% distance between the left edge of a window and the left edge of the next
% window; and corrthresh the correlation detection threshold.
%
% Advisory: let winsize = 1000 and jump = 500, corrthresh should be no
% greater than 0.075 for optimal results.
%
% Subordinate functions: <findautomin>

fast = 0;

P = zeros( ((1 + ceil(size(s,1) / winsize)) * winsize)  , 1);
length = ceil(size(s,1) / winsize) * winsize;
P(size(P,1) : size(P,1) + 1200) = 0;
W = 0;
Tprev = -1;
windownumber = 1;

s(size(s, 1) : (ceil(size(s,1) / winsize) * winsize ) + 1200) = 0;

for i = 1 : jump : length %(ceil(size(P,1) / winsize) * winsize)
    
    if(Tprev > 40)
        % This executes if we know the pitch from the previous frame and
        % assumes that we can find the new pitch reasonably close to the
        % old one.
        
        shift = findautomin(s, i, Tprev - 20, Tprev + 20, corrthresh);
        Tprev = shift;
        fast = fast + 1;
        
    else
        % If we don't have an "old pitch", find one for this frame
        shift = findautomin(s, i, 10, 800, corrthresh);
        Tprev = shift;
    end
    
    % Calculate the frequency
    frequency = 44100 / shift;
        
    if(shift > 70)
        
        W(windownumber, 1) = frequency;
        
        % If the entire frame was zero...
        if (sum(P(i : i + winsize)) == 0)
            P(i : i + winsize) = frequency;
        else
            
            % Otherwise, take an average of what was there before
            
            for j = i : i + winsize
                
                if( P(j, 1) == 0)
                    P(j, 1) = frequency;
                else
                    P(j, 1) = (P(j, 1) + frequency) / 2;
                end
            end
            
        end
        
    else
        W(windownumber, 1) = 0;
    end

    windownumber = windownumber + 1;
    
end

%W = W(1: ceil( (size(s, 1) - winsize) / jump) + 1);
sprintf('FAST-Autocorrelated windows: %d', fast)
sprintf('Fraction of total: %g', fast / windownumber)