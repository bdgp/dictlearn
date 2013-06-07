function outputChar = blankKill(inputChar)

outputChar = inputChar;
for i = 1:length(inputChar)
    charTemp = inputChar{i};
    for j = 1:length(charTemp)
        if length(strmatch(charTemp(j),' '))==0
            outputChar{i} = charTemp(j:end);
            break;
        end
    end
end