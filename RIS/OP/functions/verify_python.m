function [status] = verify_python()

    % Verifique se o Python já está configurado
    current_py_config = pyversion;
    if ~isempty(current_py_config)
        disp('Python is already configure in MATLAB.');
        status = 0;
        return; % Saia da função ou script
    end
    
    % Verifique o sistema operacional
    if ispc
        % Windows
        command = 'where python3.9';
    else
        % Linux
        command = 'which python3';
    end
    
    % Execute o comando
    [status, output] = system(command);
    
    % Verifique se a chamada do sistema foi bem sucedida
    if status == 0
        % O caminho do Python está na variável 'output'
        python_path = strtrim(output);
        disp(['Python Path: ' python_path]);
        disp('Setting up python...');
        pyversion([python_path]);
        pyenv
        disp('Python ready');
    else
        error('Error to search python path. Exiting...');
        % error('Exiting...');
        return; % Saia da função ou script
    end
    
end
    