function extract_event_by_channel_Callback(hObject, ~)

% import package contents

clear_existing_events = 1;

handles = guidata(hObject);


[currentVal, currentName, noOfEntries, str] = kVIS_dataSetListState(hObject)


fds = kVIS_getCurrentFds(hObject);

% kVIS_updateEventList(hObject, fds.eventList, false);

for ii = 1:noOfEntries

    % Extract channel
    fprintf('Generating events based on RCIN changes.\n');
    times   = kVIS_fdsGetChannel(fds, 'RCIN','Time');
    channel = kVIS_fdsGetChannel(fds, 'RCIN','C6');

    in = [];
    out = [];

    % Loop through and capture times when channel was above 1200
    idx = find(channel > 1700);

    in(1) = idx(1);

    for jj = 2:numel(idx)
        if idx(jj) - idx(jj-1) > 1
            % We had a toggle of the channel
            out(end+1) = jj-1;
            in(end+1) = jj;
        end
    end

    % Add the last out if we finished during a manourve
    if numel(in) ~= numel(out)
        out(end+1) = idx(end);
    end

    % Get the times and fill eList
    eventNumber = 0;
    for jj = 1:numel(in)

            % Fill out eList
            eventNumber = eventNumber+1;
            eList(eventNumber).type = 'Ch6 High';
            eList(eventNumber).start= times(in(jj))-2;
            eList(eventNumber).end  = times(out(jj));
            eList(eventNumber).description = 'Ch6 High';
            eList(eventNumber).plotDef='';

    end

    % Concatenate the eLists
    if (clear_existing_events)
        fds.eventList = eList;
    else
        fds.eventList = [fds.eventList,eList];
    end

    % Update data sets
    kVIS_updateEventList(hObject, fds.eventList, false);
 
end

% All done
return

end
