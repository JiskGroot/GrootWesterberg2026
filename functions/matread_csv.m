function data = matread_csv(filename)

if isstring(filename)
    filename = char(filename);
end

data = readtable(filename);
data.Properties.VariableNames = readcell([filename(1:end-4) '_header.csv']);

end