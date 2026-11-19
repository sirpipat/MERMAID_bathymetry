% RAYPATH  Container of a seismic ray path in layered media
%
%   Each entry contains
%       name      cell array of phase in each leg of the path
%                 For example, '2upS-end' is an upgoing S-wave in the
%                 second layer from the top which is the final leg of the
%                 path, denoted by '-end' tag
%       time      total travel time of all legs
%       exptime   expected minimum travel time possible to reach end point
%       amp       amplitude
%
%   Usage
%       rp = RayPath('4upP', 1.2, 0, 1);    % intialization
%       
% Last modified by spipatprathanporn@ucsd.edu, 11/18/2026

classdef RayPath
    properties
        name
        time
        exptime
        amp
    end

    methods
        function obj = RayPath(name, time, exptime, amp)
            obj.name = {name};
            obj.time = time;
            obj.exptime = exptime;
            obj.amp = amp;
        end

        function tf = isokay(obj, tmax, ampmin)
            tf = and(obj.exptime <= tmax, abs(obj.amp) >= abs(ampmin));
        end

        function tf = isdone(obj)
            tf = contains(obj.name{end}, 'end');
        end

        function n = currlayer(obj)
            n = sscanf(obj.name{end}, '%d', 1);
        end

        function ph = currphase(obj)
            ph = obj.name{end};
        end

        function ud = currdirection(obj)
            if contains(obj.name{end}, 'up')
                ud = 'up';
            else
                ud = 'down';
            end
        end

        function k = key(obj)
            k = [obj.exptime abs(obj.amp)];
        end

        function obj = addphase(obj, name, time, ampratio)
            obj.name{length(obj.name)+1} = name;
            obj.time = obj.time + time;
            obj.amp = obj.amp * ampratio;
        end

        function obj = setexptime(obj, layers, targetlayer, tP_up, tP_down)
            curr = obj.name{end};
            currlayer = sscanf(curr, '%d', 1);
            if currlayer > targetlayer
                if contains(curr,'up')
                    currlayer = currlayer - 1;
                end
                t = sum([layers(targetlayer+1:currlayer).tP]) + tP_up;
            elseif currlayer < targetlayer
                if contains(curr,'down')
                    currlayer = currlayer + 1;
                end
                t = sum([layers(currlayer:targetlayer-1).tP]) + tP_down;
            else
                if contains(curr, 'up')
                    t = tP_down;
                else
                    t = tP_up;
                end
            end
            obj.exptime = obj.time + t;
        end

        function x = getx(obj, layers)
            curr = obj.name{end};
            currlayer = sscanf(curr, '%d', 1);
            if contains(curr, 'P')
                x = obj.amp * sin(layers(currlayer).theta_p);
            else
                if contains(curr, 'up')
                    x = obj.amp * -cos(layers(currlayer).theta_s);
                else
                    x = obj.amp * cos(layers(currlayer).theta_s);
                end
            end
        end

        function z = getz(obj, layers)
            curr = obj.name{end};
            currlayer = sscanf(curr, '%d', 1);
            if contains(curr, 'P')
                if contains(curr, 'up')
                    z = obj.amp * cos(layers(currlayer).theta_p);
                else
                    z = obj.amp * -cos(layers(currlayer).theta_p);
                end
            else
                z = obj.amp * -sin(layers(currlayer).theta_s);
            end
        end

        function ph = phasename(obj, istex)
            defval('istex', false)
            % initial phase
            if contains(obj.name{1}, 'P')
                if contains(obj.name{1}, 'up')
                    ph = 'p';
                else
                    ph = 'P';
                end
            else
                if contains(obj.name{1}, 'up')
                    ph = 's';
                else
                    ph = 'S';
                end
            end
            ph_curr = ph;

            % loop through legs
            for ii = 2:length(obj.name)
                % check the next phase
                if contains(obj.name{ii}, 'P')
                    if contains(obj.name{ii}, 'up')
                        ph_next = 'p';
                    else
                        ph_next = 'P';
                    end
                else
                    if contains(obj.name{ii}, 'up')
                        ph_next = 's';
                    else
                        ph_next = 'S';
                    end
                end
                % layer of the next phase
                n = sscanf(obj.name{ii}, '%d', 1);
                
                % same propagation direction
                if (isstrprop(ph_curr, 'upper') == ...
                        isstrprop(ph_next, 'upper'))
                    % check if it changes from P to S or vice-versa
                    if ~strcmp(ph_curr, ph_next)
                        % going up
                        if isstrprop(ph_next, 'lower')
                            ph = sprintf('%s%d%s', ph, n+1, ph_next);
                        % going down
                        else
                            ph = sprintf('%s%d%s', ph, n, ph_next);
                        end
                    end
                % change propagation direction
                else
                    % reflection below the top interface of the next layer
                    if isstrprop(ph_next, 'upper')
                        if istex
                            ph = sprintf('%s\\^^%d%s', ph, n, ph_next);
                        else
                            ph = sprintf('%s^%d%s', ph, n, ph_next);
                        end
                    % reflection above the bottom interface of the next layer
                    else
                        if istex
                            ph = sprintf('%sv_%d%s', ph, n+1, ph_next);
                        else
                            ph = sprintf('%sv%d%s', ph, n+1, ph_next);
                        end
                    end
                end

                % update the current phase
                ph_curr = ph_next;
            end
        end
    end
end