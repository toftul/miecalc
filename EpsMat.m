% function epsilon=EpsMat(w,Material) 
% Returns an array of complex dielectric permittivities. The input paramters 
% are photon energy w in eV, and a string variable Material. Material can
% be Material = 'Au', 'Ag', 'Si', 'Dr'
% 'Dr' is for Drude model, which is defined by the formula
% epsilon=epsinf--(wp)^2./(w.*(w+1i*gamma));

function epsilon=EpsMat(w, material)
    const = 1.23984;  % [\mu m * eV]
    wl = const./w;  % [\mu m]
    global epsinf wp gamma
    if strcmp(material, 'Si')
        data=dlmread ('refractiveindex/SiGreen.txt');
        % depends on file from refractiveindex.info
        numRows = 121;
        % first numRows elements of the first column
        wlData = data(1:numRows, 1);  % [\mum]
        nData = data(1:numRows, 2);
        kData = data(numRows+1:end, 2);
        % n,k -> epsRe, epsIm
        epsReData = nData.*nData - kData.*kData;
        epsImData = 2.*nData.*kData;
        % interpolate
        epsRe=interp1(wlData, epsReData, wl, 'linear', 1);
        epsIm=interp1(wlData, epsImData, wl, 'linear', 0);
        epsilon=complex(epsRe, epsIm);
    elseif strcmp(material, 'Ag')
        data=dlmread ('refractiveindex/AgWerner.txt');
        % depends on file from refractiveindex.info
        numRows = 150;
        % first numRows elements of the first column
        wlData = data(1:numRows, 1);  % [\mum]
        nData = data(1:numRows, 2);
        kData = data(numRows+1:end, 2);
        % n,k -> epsRe, epsIm
        epsReData = nData.*nData - kData.*kData;
        epsImData = 2.*nData.*kData;
        % interpolate
        epsRe=interp1(wlData, epsReData, wl, 'linear', 1);
        epsIm=interp1(wlData, epsImData, wl, 'linear', 0);
        epsilon=complex(epsRe, epsIm);
    elseif strcmp(material, 'Au')
        data=dlmread ('refractiveindex/AuJohnson.txt');
        % depends on file from refractiveindex.info
        numRows = 49;
        % first numRows elements of the first column
        wlData = data(1:numRows, 1);  % [\mum]
        nData = data(1:numRows, 2);
        kData = data(numRows+1:end, 2);
        % n,k -> epsRe, epsIm
        epsReData = nData.*nData - kData.*kData;
        epsImData = 2.*nData.*kData;
        % interpolate
        epsRe=interp1(wlData, epsReData, wl, 'linear', 1);
        epsIm=interp1(wlData, epsImData, wl, 'linear', 0);
        epsilon=complex(epsRe, epsIm);
    elseif strcmp(material, 'Cu')
        data=dlmread ('refractiveindex/CuJohnson.txt');
        % depends on file from refractiveindex.info
        numRows = 49;
        % first numRows elements of the first column
        wlData = data(1:numRows, 1);  % [\mum]
        nData = data(1:numRows, 2);
        kData = data(numRows+1:end, 2);
        % n,k -> epsRe, epsIm
        epsReData = nData.*nData - kData.*kData;
        epsImData = 2.*nData.*kData;
        % interpolate
        epsRe=interp1(wlData, epsReData, wl, 'linear', 1);
        epsIm=interp1(wlData, epsImData, wl, 'linear', 0);
        epsilon=complex(epsRe, epsIm);
    elseif strcmp(material, 'Al')
        data=dlmread ('refractiveindex/AlRakic.txt');
        % depends on file from refractiveindex.info
        numRows = 206;
        % first numRows elements of the first column
        wlData = data(1:numRows, 1);  % [\mum]
        nData = data(1:numRows, 2);
        kData = data(numRows+1:end, 2);
        % n,k -> epsRe, epsIm
        epsReData = nData.*nData - kData.*kData;
        epsImData = 2.*nData.*kData;
        % interpolate
        epsRe=interp1(wlData, epsReData, wl, 'linear', 1);
        epsIm=interp1(wlData, epsImData, wl, 'linear', 0);
        epsilon=complex(epsRe, epsIm);
    else  % Drude model
        omega = 2*pi*physconst('LightSpeed')./(wl.*1e-6);
        epsilon = epsinf - wp*wp / (omega.*(omega + 1j * gamma));
    end
return
