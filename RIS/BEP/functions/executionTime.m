function str_time = executionTime(execution_seconds)

    hours = fix(execution_seconds / 3600);
    minutes = fix(mod(execution_seconds, 3600) / 60);
    seconds = fix(mod(execution_seconds, 60));
    str_time = sprintf('%02dh%02dm%02ds', hours, minutes, seconds);

end