function id = classId(name)
    name = string(name);
    switch name
        case "NR",    id = 1;
        case "LTE",   id = 2;
        case "WLAN",  id = 3;
        case "RADAR", id = 4;
        otherwise,    id = 0;
    end
end