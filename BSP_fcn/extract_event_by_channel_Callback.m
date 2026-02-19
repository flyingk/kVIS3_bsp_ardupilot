function extract_event_by_channel_Callback(hObject, ~)

% import package contents

clear_existing_events = 1;

% handles = guidata(hObject);


[fds, names, ~] = kVIS_getAllFds(hObject);

prompts = {'Group','Channel','High Value', 'Event Start Offset [ s ]'};
dlgtitle = 'Extract Events';
dims = [ 1 35 ];
definput = {'RCIN','C6', '1700', '0.0'};

answer = inputdlg(prompts,dlgtitle,dims,definput);

if isempty(answer)
    fprintf('Input cancelled\n');
    return;
end

group = answer{1};
channel = answer{2};
highVal = str2double(answer{3});
startOffset = str2double(answer{4});


% Loop through the fds files
for ii = 1:numel(fds)

    % Extract channel
    fprintf('Generating events based on RCIN changes.\n');
    times   = kVIS_fdsGetChannel(fds{ii}, group,'Time');
    values = kVIS_fdsGetChannel(fds{ii}, group, channel);

    in = [];
    out = [];

    % Loop through and capture times when channel was above 1200
    idx = find(values > highVal);

    if isempty(idx)
        fprintf('\tNo events found for fds %d\n',ii);
        continue
    end

    in(1) = idx(1);

    for jj = 2:numel(idx)
        if idx(jj) - idx(jj-1) > 1
            % We had a toggle of the channel
            out(end+1) = idx(jj-1);
            in(end+1) = idx(jj);
        end
    end

    % Add the last out if we finished during a manourve
    if numel(in) ~= numel(out)
        out(end+1) = idx(end);
    end

    % Get the times and fill eList
    eventNumber = 0;
    eList = [];

    for jj = 1:numel(in)

            % Fill out eList
            eventNumber = eventNumber+1;
            eList(eventNumber).type = [channel, ' High'];
            eList(eventNumber).start= times(in(jj))+startOffset;
            eList(eventNumber).end  = times(out(jj));
            eList(eventNumber).description = [channel, ' High'];
            eList(eventNumber).plotDef='';

    end

    % Concatenate the eLists
    if (clear_existing_events)
        fds{ii}.eventList = eList;
    else
        fds{ii}.eventList = [fds{ii}.eventList,eList];
    end

    % Update data sets
    kVIS_updateDataSet(hObject, fds{ii}, names{ii});
 
end

% All done
return

end
