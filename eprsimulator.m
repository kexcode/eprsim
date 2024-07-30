classdef eprsimulator < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        EPRSimulatorUIFigure            matlab.ui.Figure
        FileMenu                        matlab.ui.container.Menu
        LoadMenu                        matlab.ui.container.Menu
        SaveMenu                        matlab.ui.container.Menu
        ExitMenu                        matlab.ui.container.Menu
        HelpMenu                        matlab.ui.container.Menu
        AboutMenu                       matlab.ui.container.Menu
        EasyspinMenu                    matlab.ui.container.Menu
        ContactAuthorMenu               matlab.ui.container.Menu
        GridLayoutMain                  matlab.ui.container.GridLayout
        GridLayoutRightCol              matlab.ui.container.GridLayout
        CWspectrumPanel                 matlab.ui.container.Panel
        GridLayoutCW                    matlab.ui.container.GridLayout
        UIAxesCW                        matlab.ui.control.UIAxes
        AbsorptionspectrumPanel         matlab.ui.container.Panel
        GridLayoutAbsorption            matlab.ui.container.GridLayout
        UIAxesAbsorption                matlab.ui.control.UIAxes
        EnergyleveldiagramPanel         matlab.ui.container.Panel
        GridLayoutLevels                matlab.ui.container.GridLayout
        UIAxesLevels                    matlab.ui.control.UIAxes
        Panel                           matlab.ui.container.Panel
        GridLayoutLeftPanel             matlab.ui.container.GridLayout
        ControlPanel                    matlab.ui.container.Panel
        GridLayout10                    matlab.ui.container.GridLayout
        GridLayout13                    matlab.ui.container.GridLayout
        ClearButton                     matlab.ui.control.Button
        DataisloadedLamp                matlab.ui.control.Lamp
        DataisloadedLampLabel           matlab.ui.control.Label
        GridLayout12                    matlab.ui.container.GridLayout
        StartButton                     matlab.ui.control.StateButton
        GridLayout11_2                  matlab.ui.container.GridLayout
        OverlayresultsSwitch            matlab.ui.control.Switch
        OverlayresultsSwitchLabel       matlab.ui.control.Label
        GridLayout11                    matlab.ui.container.GridLayout
        SingleorientationSwitch         matlab.ui.control.Switch
        SingleorientationSwitchLabel    matlab.ui.control.Label
        GridLayout7                     matlab.ui.container.GridLayout
        ExamplesPanel                   matlab.ui.container.Panel
        GridLayout9                     matlab.ui.container.GridLayout
        ExamplesDropDown                matlab.ui.control.DropDown
        FrequencyBandPanel              matlab.ui.container.Panel
        GridLayout8                     matlab.ui.container.GridLayout
        FrequencyBandDropDown           matlab.ui.control.DropDown
        MessageField                    matlab.ui.control.EditField
        ExperimentalParametersPanel     matlab.ui.container.Panel
        GridLayoutExp                   matlab.ui.container.GridLayout
        EuleranglesdegLabel             matlab.ui.control.Label
        TemperatureKLabel               matlab.ui.control.Label
        NumberofpointsLabel             matlab.ui.control.Label
        ExpCrystalOrientationEditField  matlab.ui.control.EditField
        ExpCrystalOrientationEditFieldLabel  matlab.ui.control.Label
        ExpTemperatureEditField         matlab.ui.control.EditField
        ExpTemperatureEditFieldLabel    matlab.ui.control.Label
        ExpnPointsEditField             matlab.ui.control.NumericEditField
        ExpnPointsEditFieldLabel        matlab.ui.control.Label
        FieldrangemTLabel               matlab.ui.control.Label
        mwfrequencyGHzLabel             matlab.ui.control.Label
        ExpRangeEditField               matlab.ui.control.EditField
        ExpRangeEditFieldLabel          matlab.ui.control.Label
        ExpmwFreqEditField              matlab.ui.control.NumericEditField
        ExpmwFreqEditFieldLabel         matlab.ui.control.Label
        SpinSystemPanel                 matlab.ui.container.Panel
        GridLayoutSpinSystem            matlab.ui.container.GridLayout
        SyslwEditField                  matlab.ui.control.EditField
        SyslwEditFieldLabel             matlab.ui.control.Label
        SysSEditField                   matlab.ui.control.EditField
        SysSEditFieldLabel              matlab.ui.control.Label
        GridLayoutGn                    matlab.ui.container.GridLayout
        NucGnVal                        matlab.ui.control.Label
        gnLabel                         matlab.ui.control.Label
        GridLayoutNucSpin               matlab.ui.container.GridLayout
        NucSpinVal                      matlab.ui.control.Label
        ILabel                          matlab.ui.control.Label
        LinewidthmTLabel                matlab.ui.control.Label
        ZFScm1Label                     matlab.ui.control.Label
        QuadrupoletensorMHzLabel        matlab.ui.control.Label
        HyperfinetensorMHzLabel         matlab.ui.control.Label
        SysDEditField                   matlab.ui.control.EditField
        SysDEditFieldLabel              matlab.ui.control.Label
        SysQEditField                   matlab.ui.control.EditField
        SysQEditFieldLabel              matlab.ui.control.Label
        SysAEditField                   matlab.ui.control.EditField
        SysAEditFieldLabel              matlab.ui.control.Label
        NucleargvalueLabel              matlab.ui.control.Label
        NuclearspinLabel                matlab.ui.control.Label
        ofequivalentnucleiLabel         matlab.ui.control.Label
        NucleiLabel                     matlab.ui.control.Label
        SysnEditField                   matlab.ui.control.EditField
        SysnEditFieldLabel              matlab.ui.control.Label
        SysNucsEditField                matlab.ui.control.EditField
        SysNucsEditFieldLabel           matlab.ui.control.Label
        SysgEditField                   matlab.ui.control.EditField
        SysgEditFieldLabel              matlab.ui.control.Label
        ElectronspingtensorLabel        matlab.ui.control.Label
        ElectronspinLabel               matlab.ui.control.Label
        ComputationalParametersPanel    matlab.ui.container.Panel
        GridLayoutOpt                   matlab.ui.container.GridLayout
        NumberofknotsLabel              matlab.ui.control.Label
        CalculationmethodLabel          matlab.ui.control.Label
        OptnKnotsEditField              matlab.ui.control.NumericEditField
        OptnKnotsEditFieldLabel         matlab.ui.control.Label
        OptMethodDropDown               matlab.ui.control.DropDown
        OptMethodDropDownLabel          matlab.ui.control.Label
    end

    
    properties (Access = private)
        Sys = struct('S', 0.5, 'g', 2.0023, 'Nucs', '','lw', 0.1);          % Spin System parameters
        Exp = struct('Temperature', [], 'CrystalOrientation', []);          % Experimental parameters
        Opt = struct('Method', 'perturb2', 'nKnots', 23);                   % Optimization parameters
        Data = struct('x',[],'spec',[],'par',[]);                           % Loaded experimental data (optional)
        
        ColorScheme = struct('Map', lines(10), 'CurrentColorIndex', 1);
        version = '3.2'; % version
        Bohr_magn = 0.927400968e-20; % [Erg/G]
        planck_const = 2 * pi * 1.054571726e-27; % 2pi*planck const = h [Erg*s] 
        C = 2.99792458e10; % [cm/s]
        %gamma_e = 2.0023 * app.Bohr_magn * 10 / planck_const;    % in Hz/mT
    end
    
    methods (Access = private)
        
        function update_spin_parameters(app)
            if (~isempty(app.SysSEditField.Value))
                app.Sys.S = str2num(app.SysSEditField.Value);
            elseif (isfield(app.Sys,'S'))
                app.Sys = rmfield(app.Sys, 'S');
            end
            if (~isempty(app.SysgEditField.Value))
                app.Sys.g = str2num(app.SysgEditField.Value);
            elseif (isfield(app.Sys,'G'))
                app.Sys = rmfield(app.Sys, 'G');
            end            
            if (~isempty(app.SysNucsEditField.Value))
                app.Sys.Nucs = app.SysNucsEditField.Value;
            elseif (isfield(app.Sys,'Nucs'))
                app.Sys = rmfield(app.Sys, 'Nucs');
            end
            if (~isempty(app.SysnEditField.Value))
                app.Sys.n = str2num(app.SysnEditField.Value);
            elseif (isfield(app.Sys,'n'))
                app.Sys = rmfield(app.Sys, 'n');
            end
            if (~isempty(app.SysAEditField.Value))
                app.Sys.A = str2num(app.SysAEditField.Value);
            elseif (isfield(app.Sys,'A'))
                app.Sys = rmfield(app.Sys, 'A');
            end            
            if (~isempty(app.SysQEditField.Value))
                app.Sys.Q = str2num(app.SysQEditField.Value);
            elseif (isfield(app.Sys,'Q'))
                app.Sys = rmfield(app.Sys, 'Q');
            end            
            if (~isempty(app.SysDEditField.Value))
                app.Sys.D = app.C / 1e6 * str2num(app.SysDEditField.Value);
            elseif (isfield(app.Sys,'D'))
                app.Sys = rmfield(app.Sys, 'D');
            end            
            if (~isempty(app.SyslwEditField.Value))
                app.Sys.lw = str2num(app.SyslwEditField.Value);
            elseif (isfield(app.Sys,'lw'))
                app.Sys = rmfield(app.Sys, 'lw');
            end            
        end
                
        function update_exp_parameters(app)
            if (~isempty(app.ExpmwFreqEditField.Value))
                app.Exp.mwFreq = app.ExpmwFreqEditField.Value;
            elseif (isfield(app.Exp,'mwFreq'))
                app.Exp = rmfield(app.Exp, 'mwFreq');
            end
            if (~isempty(app.ExpRangeEditField.Value))
                app.Exp.Range = str2num(app.ExpRangeEditField.Value);
            elseif (isfield(app.Exp,'Range'))
                app.Exp = rmfield(app.Exp, 'Range');
            end
            if (~isempty(app.ExpnPointsEditField.Value))
                app.Exp.nPoints = app.ExpnPointsEditField.Value;
            elseif (isfield(app.Exp,'nPoints'))
                app.Exp = rmfield(app.Exp, 'nPoints');
            end
            if (~isempty(app.ExpTemperatureEditField.Value))
                app.Exp.Temperature = str2num(app.ExpTemperatureEditField.Value);
            elseif (isfield(app.Exp,'Temperature'))
                app.Exp = rmfield(app.Exp, 'Temperature');
            end
            if (~isempty(app.ExpCrystalOrientationEditField.Value) && strcmp(app.SingleorientationSwitch.Value, 'On'))
                app.Exp.CrystalOrientation = str2num(app.ExpCrystalOrientationEditField.Value);    
            elseif (isfield(app.Exp,'CrystalOrientation'))
                app.Exp = rmfield(app.Exp, 'CrystalOrientation');
            end
        end
        
        function update_opt_parameters(app)
            app.Opt.Method = app.OptMethodDropDown.Value;
            if (~isempty(app.OptnKnotsEditField.Value))
                app.Opt.nKnots = app.OptnKnotsEditField.Value;
            elseif (isfield(app.Opt,'nKnots'))
                app.Opt = rmfield(app.Opt, 'nKnots');
            end
        end
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Value changed function: SysgEditField
        function SysgEditFieldValueChanged(app, event)
            app.Sys.g = str2num(app.SysgEditField.Value);
        end

        % Value changed function: SysNucsEditField
        function SysNucsEditFieldValueChanged(app, event)
            value = app.SysNucsEditField.Value;
            try
                nucdata(value);
                app.NucSpinVal.Text = num2str(nucspin(value));
                app.NucGnVal.Text = num2str(nucgval(value));
                app.MessageField.Value = '';
            catch ME
                if (~isempty(ME.message))
                    if (contains(ME.message, 'Unknown element'))
                        app.MessageField.Value = ['Error: ' ME.message];
                    end
                    app.NucSpinVal.Text = '';
                    app.NucGnVal.Text = '';
                end
            end
        end

        % Value changed function: SysnEditField
        function SysnEditFieldValueChanged(app, event)
            if (strcmp(app.SysNucsEditField.Value,''))
                app.SysnEditField.Value = '';
            else
                app.Sys.n = round(str2num(app.SysnEditField.Value));
            end
        end

        % Menu selected function: LoadMenu
        function LoadMenuSelected(app, event)
            % select file and load EPR data
            warning('off'); 
            try
                [app.Data.file, app.Data.path] = uigetfile('*.dsc');
                [app.Data.x, app.Data.spec, app.Data.par] = eprload(strcat(app.Data.path,app.Data.file)); % loading selected data file | only .dsc are supported for now
                app.MessageField.FontColor = 'k';            
                app.MessageField.Value = "Loading...";
                
                app.Data.x = app.Data.x/10; % converting G into mT
                app.Data.spec = app.Data.spec/max(app.Data.spec); % scaling ordinate to maxval == 1
                if (isfield(app.Data.par, 'MWFQ'))
                    % change band selector depending on mwFreq
                    freq = app.Data.par.MWFQ/1e9
                    if (freq > 2 && freq < 4)
                        app.FrequencyBandDropDown.Value = 'S-band';
                    elseif (freq > 30 && freq < 50)
                        app.FrequencyBandDropDown.Value = 'Q-band';
                    elseif (freq > 75 && freq < 110)
                        app.FrequencyBandDropDown.Value = 'W-band';
                    elseif (freq > 230 && freq < 250)
                        app.FrequencyBandDropDown.Value = '240 GHz';
                    else
                        app.FrequencyBandDropDown.Value = 'X-band';
                    end
                    app.ExpmwFreqEditField.Value = app.Data.par.MWFQ/1e9; % loading frequency if exists    
                end
                app.Data.par.Range = [app.Data.par.XMIN (app.Data.par.XMIN+app.Data.par.XWID)]/10; % loading field range
                app.ExpRangeEditField.Value = strcat('[',num2str(app.Data.par.Range(1))," ",num2str(app.Data.par.Range(2)),']'); 
      
                % is pulse or cw experiment?
                if (isfield(app.Data.par, 'EXPT'))
                    if (strcmp(app.Data.par.EXPT,'PLS'))
                        if (isreal(app.Data.spec)) % check if not complex
                            plot(app.UIAxesAbsorption, app.Data.x, app.Data.spec,'k')
                        else
                            plot(app.UIAxesAbsorption, app.Data.x, real(app.Data.spec),'k',app.Data.x, imag(app.Data.spec),'m')
                        end
                        app.UIAxesAbsorption.YLim = [-0.1 1.2];
                    end
                else % otherwise / not pulse
                    plot(app.UIAxesCW, app.Data.x, app.Data.spec,'k')
                    app.UIAxesCW.YLim = 1.2 * [min(real(app.Data.spec)) max(real(app.Data.spec))];    
                end
                app.UIAxesCW.XLim = app.Data.par.Range;
                app.UIAxesAbsorption.XLim = app.Data.par.Range;
                app.UIAxesLevels.XLim = app.Data.par.Range;                
                
                % lamp on
                app.DataisloadedLamp.Color = 'g';
                % status message
                app.MessageField.Value = "Loading completed";
            catch ME
                app.MessageField.Value = ME.message;
            end
            app.EPRSimulatorUIFigure.Visible = 'off';
            app.EPRSimulatorUIFigure.Visible = 'on';
            warning('on');
        end

        % Value changed function: SysSEditField
        function SysSEditFieldValueChanged(app, event)
            app.Sys.S = str2num(app.SysSEditField.Value);   
        end

        % Menu selected function: SaveMenu
        function SaveMenuSelected(app, event)
            update_spin_parameters(app);
            update_exp_parameters(app);
            update_opt_parameters(app);
            [file, path] = uiputfile('*.mat','File Selection','parameters.mat');
            Sys = app.Sys; 
            save(strcat(path, file), 'Sys');
            Exp = app.Exp; 
            save(strcat(path, file), 'Exp', '-append');
            Opt = app.Opt; 
            save(strcat(path, file), 'Opt', '-append');
            clear Sys Exp Opt;
            app.MessageField.Value = "Parameters saved";
        end

        % Value changed function: SysAEditField
        function SysAEditFieldValueChanged(app, event)
            app.Sys.A = str2num(app.SysAEditField.Value);
        end

        % Value changed function: SysQEditField
        function SysQEditFieldValueChanged(app, event)
            app.Sys.Q = str2num(app.SysQEditField.Value);
        end

        % Value changed function: SysDEditField
        function SysDEditFieldValueChanged(app, event)
            app.Sys.D = app.C / 1e6 * str2num(app.SysDEditField.Value);
        end

        % Value changed function: SyslwEditField
        function SyslwEditFieldValueChanged(app, event)
            app.Sys.lw = str2num(app.SyslwEditField.Value);
        end

        % Value changed function: ExpmwFreqEditField
        function ExpmwFreqEditFieldValueChanged(app, event)
            app.Exp.mwFreq = app.ExpmwFreqEditField.Value;
        end

        % Value changed function: ExpRangeEditField
        function ExpRangeEditFieldValueChanged(app, event)
            app.Exp.Range = str2num(app.ExpRangeEditField.Value);
        end

        % Value changed function: ExpnPointsEditField
        function ExpnPointsEditFieldValueChanged(app, event)
            app.Exp.nPoints = app.ExpnPointsEditField.Value;
        end

        % Value changed function: ExpTemperatureEditField
        function ExpTemperatureEditFieldValueChanged(app, event)
            app.Exp.Temperature = str2num(app.ExpTemperatureEditField.Value);
        end

        % Value changed function: ExpCrystalOrientationEditField
        function ExpCrystalOrientationEditFieldValueChanged(app, event)
            app.Exp.CrystalOrientation = pi / 180 * str2num(app.ExpCrystalOrientationEditField.Value);
        end

        % Value changed function: OptMethodDropDown
        function OptMethodDropDownValueChanged(app, event)
            app.Opt.Method = app.OptMethodDropDown.Value;
        end

        % Value changed function: OptnKnotsEditField
        function OptnKnotsEditFieldValueChanged(app, event)
            app.Opt.nKnots = round(app.OptnKnotsEditField.Value);
        end

        % Value changed function: SingleorientationSwitch
        function SingleorientationSwitchValueChanged(app, event)
            oriSwitch = app.SingleorientationSwitch.Value;
            if (strcmp(oriSwitch,'On'))
                app.ExpCrystalOrientationEditField.Enable = true;
                if (isempty(app.ExpCrystalOrientationEditField.Value))
                    app.ExpCrystalOrientationEditField.Value = '[0 0 0]';
                end
            else
                app.ExpCrystalOrientationEditField.Enable = false;
                app.ExpCrystalOrientationEditField.Value = '[0 0 0]';
            end
        end

        % Value changed function: StartButton
        function StartButtonValueChanged(app, event)
            value = app.StartButton.Value;
            if (value == true)
                app.StartButton.Text = "Stop";
                app.MessageField.Value = "Busy...";
                update_spin_parameters(app);
                update_exp_parameters(app);
                update_opt_parameters(app);
                
                % calculating fieldStep
                fieldStep = (app.Exp.Range(2) - app.Exp.Range(1)) / (app.Exp.nPoints - 1);
                if(fieldStep > app.Sys.lw / 3)
                    app.MessageField.Value = "Error: not enough points to resolve the resonance line";
                    app.StartButton.Value = false;
                    app.StartButton.Text = "Start";
                    return;
                end
                
                % Energy Level Diagram
                B = linspace(app.Exp.Range(1), app.Exp.Range(2), app.Exp.nPoints)';
                app.SingleorientationSwitch.Value
                if (strcmp(app.SingleorientationSwitch.Value,'On'))
                    Ori = [app.Exp.CrystalOrientation(1) app.Exp.CrystalOrientation(2)]*pi/180;
                else
                    Ori = [0 0];
                end

                % calculate transitions
                Param = app.Exp;
                Param.CrystalOrientation = [Ori 0];
                
                ELD_Sys = app.Sys;
                
%---------------equivalent nuclei case-------------------------------------
                if (isfield(app.Sys, 'n') && app.Sys.n > 1)
                    for i = 1:app.Sys.n - 1
                        ELD_Sys.Nucs = strcat(ELD_Sys.Nucs,',',app.Sys.Nucs);
                        ELD_Sys.A = [ELD_Sys.A; app.Sys.A];
                        if (isfield(app.Sys, 'Q') && app.Sys.Q ~= 0)
                            ELD_Sys.Q = [ELD_Sys.Q; app.Sys.Q];    
                        end
                    end
                    ELD_Sys = rmfield(ELD_Sys,'n');
%                     ELD_Sys
%                     Param
                end
%---------------calculating E levels and resonant positions----------------                 
                try
                    Elev = levels(ELD_Sys, Ori, B);
                    [Pos, ampl, ~, trans] = resfields(ELD_Sys, Param);
                catch ME
                    app.MessageField.Value = ME.message;
                    app.StartButton.Value = false;
                    app.StartButton.Text = "Start";
                    return;
                end 

                Elev_cm = Elev * 1e6 / app.C;
%                 Nlev = size(Elev,2);
                
                
                % check if need to overlap spectra
                if (strcmp(app.OverlayresultsSwitch.Value,'On'))
                    hold(app.UIAxesLevels, 'on');
                    app.ColorScheme.CurrentColorIndex = mod(app.ColorScheme.CurrentColorIndex + 1, 10) + 1; % numbers in range [1, 10]
                    trcolor = app.ColorScheme.Map(app.ColorScheme.CurrentColorIndex, :);
                else
                    hold(app.UIAxesLevels, 'off');
                    trcolor='r';
                end                
                
                % plot Energy Levels vs field
                plot(app.UIAxesLevels, B, Elev_cm, 'b')
                hold(app.UIAxesLevels, 'on');
                xlim(app.UIAxesLevels, app.Exp.Range);
                
                if (size(Pos) > 0)
                    
%                     rfield = [];
%                     amp = [];
                    amplitudeTreshold = 0.1;
                    maxAmplitude = max(ampl);
                    ampl = ampl / maxAmplitude;
                    for k = 1:size(Pos,1)
                        resField = Pos(k);
%                         rfield = [rfield; resField];
%                         amp = [amp; ampl(k)];
                        [~, index] = min(abs(B - Pos(k)));
                        E1 = Elev_cm(index, trans(k, 1));                        
                        E2 = Elev_cm(index, trans(k, 2));
                        xx=[resField resField];
                        yy=[E1 E2];
                        
                        % plot only allowed transitions                        
                        if (ampl(k) > amplitudeTreshold)
                            plot(app.UIAxesLevels, xx, yy,'Color', trcolor)
                        end
                    end
%                     rfield
%                     amp
                end

% -------------------------------------------------------------------------
                % Absorption spectrum
                app.Exp.Harmonic = 0;
                try
                    [x, spec] = pepper(ELD_Sys, app.Exp, app.Opt);
                catch ME
                    if (~isempty(ME.message))
                        app.MessageField.Value = ME.message;
                        app.StartButton.Value = false;
                        app.StartButton.Text = "Start";
                        return;
                    end
                end
                spec = spec/max(spec);
                
                % check if need to overlap spectra
                if (strcmp(app.OverlayresultsSwitch.Value,'On'))
                    hold(app.UIAxesAbsorption, 'on');
                    plot(app.UIAxesAbsorption, x, spec, 'Color', trcolor);
                else
                    hold(app.UIAxesAbsorption, 'off');
                    plot(app.UIAxesAbsorption, x, spec, 'r');
                end
                
                
                xlim(app.UIAxesAbsorption, app.Exp.Range);
%                 ylim(app.UIAxesAbsorption, [-0.2 1.1*max(spec)])
                
                % if data exists - add to the plot
                if (isfield(app.Data,'file'))
                    if (isfield(app.Data.par, 'EXPT') && strcmp(app.Data.par.EXPT,'PLS'))
                        hold(app.UIAxesAbsorption, 'on');
                        if isreal(app.Data.spec)
                            plot(app.UIAxesAbsorption, app.Data.x, app.Data.spec,'k')
                        else
                            plot(app.UIAxesAbsorption, app.Data.x, real(app.Data.spec),'k',app.Data.x, imag(app.Data.spec),'m')    
                        end
                    end
                end
                
                % CW spectrum
                app.Exp.Harmonic = 1;
                try
                    [x, spec] = pepper(ELD_Sys, app.Exp, app.Opt);
                catch ME
                    if (~isempty(ME.message))
                        app.MessageField.Value = ME.message;
                        app.StartButton.Value = false;
                        app.StartButton.Text = "Start";
                        return;
                    end
                end
                spec = spec/max(spec);
                
                % check if need to overlap spectra
                if (strcmp(app.OverlayresultsSwitch.Value,'On'))
                    hold(app.UIAxesCW, 'on');
                    plot(app.UIAxesCW, x, spec, 'Color', trcolor);
                else
                    hold(app.UIAxesCW, 'off');
                    plot(app.UIAxesCW, x, spec, 'r');
                end
                xlim(app.UIAxesCW, app.Exp.Range);
                ylim(app.UIAxesCW, 1.2*[min(spec) max(spec)]);
                               
                % if data exists - add to the plot
                if (isfield(app.Data,'file'))
                    if (isfield(app.Data.par, 'EXPT'))
                        if (~strcmp(app.Data.par.EXPT,'PLS'))
                            hold(app.UIAxesCW, 'on');
                            plot(app.UIAxesCW, app.Data.x, app.Data.spec,'k');
                        end
                    else
                        hold(app.UIAxesCW, 'on');
                        plot(app.UIAxesCW, app.Data.x, app.Data.spec,'k');
                    end
                end
                
                app.MessageField.Value = "Done!";
                app.StartButton.Value = false;    
            end
            app.StartButton.Text = "Start";
        end

        % Value changed function: ExamplesDropDown
        function ExamplesDropDownValueChanged(app, event)
            value = app.ExamplesDropDown.Value;
            switch value
                case 'Free electron'
                    app.SysSEditField.Value = '1/2';
                    app.SysgEditField.Value = '2.0023';
                    app.SysNucsEditField.Value = '';
                    app.SysnEditField.Value = '';
                    app.SysAEditField.Value = '';
                    app.SysQEditField.Value = '';
                    app.SysDEditField.Value = '';
                    app.SyslwEditField.Value = '0.3';
                    SysNucsEditFieldValueChanged(app);
                    if (~isfield(app.Data.par, 'MWFQ'))
                        app.FrequencyBandDropDown.Value = 'X-band';
                        app.ExpmwFreqEditField.Value = 9.4;
                    end
                    if (~isfield(app.Data.par, 'Range'))
                        app.ExpRangeEditField.Value = '[330 340]';
                    end
                    app.ExpnPointsEditField.Value = 501;
                    app.ExpTemperatureEditField.Value = '';
                    app.ExpCrystalOrientationEditField.Value = '[0 0 0]';
                    app.OptMethodDropDown.Value = 'perturb2'; 
                    app.OptnKnotsEditField.Value = 23;
                    app.MessageField.Value = 'Free electron example selected';
                case '1 Proton'
                    app.SysSEditField.Value = '1/2';
                    app.SysgEditField.Value = '2.0023';
                    app.SysNucsEditField.Value = '1H';
                    app.SysnEditField.Value = '';
                    app.SysAEditField.Value = '20';
                    app.SysQEditField.Value = '';
                    app.SysDEditField.Value = '';
                    app.SyslwEditField.Value = '0.2';
                    SysNucsEditFieldValueChanged(app);
                    if (~isfield(app.Data.par, 'MWFQ'))
                        app.FrequencyBandDropDown.Value = 'X-band';
                        app.ExpmwFreqEditField.Value = 9.4;
                    end
                    if (~isfield(app.Data.par, 'Range'))
                        app.ExpRangeEditField.Value = '[330 340]';
                    end
                    app.ExpnPointsEditField.Value = 501;
                    app.ExpTemperatureEditField.Value = '';
                    app.ExpCrystalOrientationEditField.Value = '[0 0 0]';
                    app.OptMethodDropDown.Value = 'perturb2'; 
                    app.OptnKnotsEditField.Value = 23;
                    app.MessageField.Value = '1 Proton example selected';
                case '2 Protons'
                    app.SysSEditField.Value = '1/2';
                    app.SysgEditField.Value = '2.0023';
                    app.SysNucsEditField.Value = '1H';
                    app.SysnEditField.Value = '2';
                    app.SysAEditField.Value = '20';
                    app.SysQEditField.Value = '';
                    app.SysDEditField.Value = '';
                    app.SyslwEditField.Value = '0.2';
                    SysNucsEditFieldValueChanged(app);
                    if (~isfield(app.Data.par, 'MWFQ'))
                        app.FrequencyBandDropDown.Value = 'X-band';
                        app.ExpmwFreqEditField.Value = 9.4;
                    end
                    if (~isfield(app.Data.par, 'Range'))
                        app.ExpRangeEditField.Value = '[330 340]';
                    end
                    app.ExpnPointsEditField.Value = 501;
                    app.ExpTemperatureEditField.Value = '';
                    app.ExpCrystalOrientationEditField.Value = '[0 0 0]';
                    app.OptMethodDropDown.Value = 'perturb2'; 
                    app.OptnKnotsEditField.Value = 23;
                    app.MessageField.Value = '2 Protons example selected';
                case 'Nitroxide radical'
                    app.SysSEditField.Value = '1/2';
                    app.SysgEditField.Value = '[2.0083 2.0061 2.0022]';
                    app.SysNucsEditField.Value = '14N';
                    app.SysnEditField.Value = '';
                    app.SysAEditField.Value = '[11.2 11.2 92.4]';
                    app.SysQEditField.Value = '';
                    app.SysDEditField.Value = '';
                    app.SyslwEditField.Value = '0.5';
                    SysNucsEditFieldValueChanged(app);
                    if (~isfield(app.Data.par, 'MWFQ'))
                        app.FrequencyBandDropDown.Value = 'X-band';
                        app.ExpmwFreqEditField.Value = 9.4;
                    end
                    if (~isfield(app.Data.par, 'Range'))
                        app.ExpRangeEditField.Value = '[328 342]';
                    end
                    app.ExpnPointsEditField.Value = 501;
                    app.ExpTemperatureEditField.Value = '';
                    app.ExpCrystalOrientationEditField.Value = '[0 0 0]';
                    app.OptMethodDropDown.Value = 'perturb2'; 
                    app.OptnKnotsEditField.Value = 23;
                    app.MessageField.Value = 'Nitroxide radical example selected';
                case 'Mn(III) ion'
                    app.SysSEditField.Value = '2';
                    app.SysgEditField.Value = '[2.000116 2.0007677 2.0073025]';
                    app.SysNucsEditField.Value = '55Mn';
                    app.SysnEditField.Value = '';
                    app.SysAEditField.Value = '290';
                    app.SysQEditField.Value = '[0 0]';
                    app.SysDEditField.Value = '[-3.968 0.0]';
                    app.SyslwEditField.Value = '20';
                    SysNucsEditFieldValueChanged(app);
                    if (~isfield(app.Data.par, 'MWFQ'))
                        app.FrequencyBandDropDown.Value = 'X-band';
                        app.ExpmwFreqEditField.Value = 9.4;
                    end
                    if (~isfield(app.Data.par, 'Range'))
                        app.ExpRangeEditField.Value = '[0 1000]';
                    end
                    app.ExpnPointsEditField.Value = 201;
                    app.ExpTemperatureEditField.Value = '';
                    app.ExpCrystalOrientationEditField.Value = '[0 0 0]';
                    app.OptMethodDropDown.Value = 'matrix'; 
                    app.OptnKnotsEditField.Value = 13;
                    app.MessageField.Value = 'Mn(III) ion example selected';
                case 'Fe(III) ion'
                    app.SysSEditField.Value = '2.5';
                    app.SysgEditField.Value = '2';
                    app.SysNucsEditField.Value = '';
                    app.SysnEditField.Value = '';
                    app.SysAEditField.Value = '';
                    app.SysQEditField.Value = '';
                    app.SysDEditField.Value = '[5 0]';
                    app.SyslwEditField.Value = '20';
                    SysNucsEditFieldValueChanged(app);
                    if (~isfield(app.Data.par, 'MWFQ'))
                        app.FrequencyBandDropDown.Value = 'X-band';
                        app.ExpmwFreqEditField.Value = 9.4;
                    end
                    if (~isfield(app.Data.par, 'Range'))
                        app.ExpRangeEditField.Value = '[0 400]';
                    end
                    app.ExpnPointsEditField.Value = 501;
                    app.ExpTemperatureEditField.Value = '';
                    app.ExpCrystalOrientationEditField.Value = '[0 0 0]';
                    app.OptMethodDropDown.Value = 'matrix'; 
                    app.OptnKnotsEditField.Value = 17;
                    app.MessageField.Value = 'Fe(III) ion example selected';
                otherwise % Spin Triplet
                    app.SysSEditField.Value = '1';
                    app.SysgEditField.Value = '2';
                    app.SysNucsEditField.Value = '';
                    app.SysnEditField.Value = '';
                    app.SysAEditField.Value = '';
                    app.SysQEditField.Value = '';
                    app.SysDEditField.Value = '[0.15 0.025]';
                    app.SyslwEditField.Value = '10';
                    SysNucsEditFieldValueChanged(app);
                    if (~isfield(app.Data.par, 'MWFQ'))
                        app.FrequencyBandDropDown.Value = 'X-band';
                        app.ExpmwFreqEditField.Value = 9.3;
                    end
                    if (~isfield(app.Data.par, 'Range'))
                        app.ExpRangeEditField.Value = '[0 500]';
                    end
                    app.ExpnPointsEditField.Value = 501;
                    app.ExpTemperatureEditField.Value = '';
                    app.ExpCrystalOrientationEditField.Value = '[0 0 0]';
                    app.OptMethodDropDown.Value = 'matrix'; 
                    app.OptnKnotsEditField.Value = 17;
                    app.MessageField.Value = 'Spin Triplet example selected';
            end
        end

        % Value changed function: FrequencyBandDropDown
        function FrequencyBandDropDownValueChanged(app, event)
            value = app.FrequencyBandDropDown.Value;
            switch value
                case 'S-band'
                    app.ExpmwFreqEditField.Value = 3.4;
                    if (~isfield(app.Data.par, 'Range'))
                        app.ExpRangeEditField.Value = '[114 128]';
                    end
                    app.MessageField.Value = 'Frequency band changed to S-band';
                case 'X-band'
                    app.ExpmwFreqEditField.Value = 9.4;
                    if (~isfield(app.Data.par, 'Range'))
                        app.ExpRangeEditField.Value = '[330 340]';
                    end
                    app.MessageField.Value = 'Frequency band changed to X-band';
                case 'Q-band'
                    app.ExpmwFreqEditField.Value = 34;
                    if (~isfield(app.Data.par, 'Range'))
                        app.ExpRangeEditField.Value = '[1205 1220]';
                    end
                    app.MessageField.Value = 'Frequency band changed to Q-band';
                case 'W-band'
                    app.ExpmwFreqEditField.Value = 94;
                    if (~isfield(app.Data.par, 'Range'))
                        app.ExpRangeEditField.Value = '[3340 3362]';
                    end
                    app.MessageField.Value = 'Frequency band changed to W-band';
                otherwise
                    app.ExpmwFreqEditField.Value = 240;
                    if (~isfield(app.Data.par, 'Range'))
                        app.ExpRangeEditField.Value = '[8530 8575]';
                    end
                    app.MessageField.Value = 'Frequency band changed to 240 GHz';
            end
        end

        % Value changed function: OverlayresultsSwitch
        function OverlayresultsSwitchValueChanged(app, event)
            value = app.OverlayresultsSwitch.Value;
            if (strcmp(value,'On'))
                hold(app.UIAxesLevels, 'on');
                hold(app.UIAxesAbsorption, 'on');
                hold(app.UIAxesCW, 'on');
                app.MessageField.Value = 'Overlaying plots on';
            else
                hold(app.UIAxesLevels, 'off');
                hold(app.UIAxesAbsorption, 'off');
                hold(app.UIAxesCW, 'off');
                app.MessageField.Value = 'Overlaying plots off';
            end
        end

        % Button pushed function: ClearButton
        function ClearButtonPushed(app, event)
            app.Data = struct('');
            app.DataisloadedLamp.Color = [0.5 0.5 0.5];
            cla(app.UIAxesLevels);
            cla(app.UIAxesAbsorption);
            cla(app.UIAxesCW);
            app.MessageField.Value = 'The data is cleared from memory';
        end

        % Menu selected function: EasyspinMenu
        function EasyspinMenuSelected(app, event)
            web('https://www.easyspin.org/');
        end

        % Menu selected function: AboutMenu
        function AboutMenuSelected(app, event)
            CreateStruct.Interpreter = 'tex';
            CreateStruct.WindowStyle = 'modal';          
            message = sprintf('EPR Simulator %s\n\nAuthor: Leonid Rapatskiy\n\nThis program displays the energy level diagram and the EPR spectrum of a given spin system. All calculations are performed using the EasySpin package (easyspin.org)', app.version);
            h = msgbox(message, 'About EPR Simulator', CreateStruct);
            set(h, 'Position', [200 500 300 120])
        end

        % Menu selected function: ContactAuthorMenu
        function ContactAuthorMenuSelected(app, event)
            web('mailto:leonid.rapatskiy@gmail.com');
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create EPRSimulatorUIFigure and hide until all components are created
            app.EPRSimulatorUIFigure = uifigure('Visible', 'off');
            app.EPRSimulatorUIFigure.Color = [1 1 1];
            app.EPRSimulatorUIFigure.Position = [100 100 1200 725];
            app.EPRSimulatorUIFigure.Name = 'EPR Simulator';

            % Create FileMenu
            app.FileMenu = uimenu(app.EPRSimulatorUIFigure);
            app.FileMenu.Text = 'File';

            % Create LoadMenu
            app.LoadMenu = uimenu(app.FileMenu);
            app.LoadMenu.MenuSelectedFcn = createCallbackFcn(app, @LoadMenuSelected, true);
            app.LoadMenu.Text = 'Load';

            % Create SaveMenu
            app.SaveMenu = uimenu(app.FileMenu);
            app.SaveMenu.MenuSelectedFcn = createCallbackFcn(app, @SaveMenuSelected, true);
            app.SaveMenu.Text = 'Save';

            % Create ExitMenu
            app.ExitMenu = uimenu(app.FileMenu);
            app.ExitMenu.Text = 'Exit';

            % Create HelpMenu
            app.HelpMenu = uimenu(app.EPRSimulatorUIFigure);
            app.HelpMenu.Text = 'Help';

            % Create AboutMenu
            app.AboutMenu = uimenu(app.HelpMenu);
            app.AboutMenu.MenuSelectedFcn = createCallbackFcn(app, @AboutMenuSelected, true);
            app.AboutMenu.Text = 'About';

            % Create EasyspinMenu
            app.EasyspinMenu = uimenu(app.HelpMenu);
            app.EasyspinMenu.MenuSelectedFcn = createCallbackFcn(app, @EasyspinMenuSelected, true);
            app.EasyspinMenu.Text = 'Easyspin';

            % Create ContactAuthorMenu
            app.ContactAuthorMenu = uimenu(app.HelpMenu);
            app.ContactAuthorMenu.MenuSelectedFcn = createCallbackFcn(app, @ContactAuthorMenuSelected, true);
            app.ContactAuthorMenu.Text = 'Contact Author';

            % Create GridLayoutMain
            app.GridLayoutMain = uigridlayout(app.EPRSimulatorUIFigure);
            app.GridLayoutMain.ColumnWidth = {'1.3x', '1x', '1x'};
            app.GridLayoutMain.RowHeight = {'1x'};

            % Create Panel
            app.Panel = uipanel(app.GridLayoutMain);
            app.Panel.BorderType = 'none';
            app.Panel.BackgroundColor = [1 1 1];
            app.Panel.HandleVisibility = 'off';
            app.Panel.Layout.Row = 1;
            app.Panel.Layout.Column = 1;

            % Create GridLayoutLeftPanel
            app.GridLayoutLeftPanel = uigridlayout(app.Panel);
            app.GridLayoutLeftPanel.ColumnWidth = {'1x'};
            app.GridLayoutLeftPanel.RowHeight = {'11.1x', '6.1x', '3.1x', '2.1x', '3.1x', '1x'};
            app.GridLayoutLeftPanel.Padding = [0 0 0 0];

            % Create ComputationalParametersPanel
            app.ComputationalParametersPanel = uipanel(app.GridLayoutLeftPanel);
            app.ComputationalParametersPanel.ForegroundColor = [0.502 0.502 0.502];
            app.ComputationalParametersPanel.Title = 'Computational Parameters';
            app.ComputationalParametersPanel.BackgroundColor = [1 1 1];
            app.ComputationalParametersPanel.Layout.Row = 3;
            app.ComputationalParametersPanel.Layout.Column = 1;
            app.ComputationalParametersPanel.FontWeight = 'bold';

            % Create GridLayoutOpt
            app.GridLayoutOpt = uigridlayout(app.ComputationalParametersPanel);
            app.GridLayoutOpt.ColumnWidth = {'2.5x', '3x', '3x'};
            app.GridLayoutOpt.RowSpacing = 5;
            app.GridLayoutOpt.Padding = [2 2 2 2];

            % Create OptMethodDropDownLabel
            app.OptMethodDropDownLabel = uilabel(app.GridLayoutOpt);
            app.OptMethodDropDownLabel.HorizontalAlignment = 'center';
            app.OptMethodDropDownLabel.Layout.Row = 1;
            app.OptMethodDropDownLabel.Layout.Column = 1;
            app.OptMethodDropDownLabel.Text = 'Opt.Method';

            % Create OptMethodDropDown
            app.OptMethodDropDown = uidropdown(app.GridLayoutOpt);
            app.OptMethodDropDown.Items = {'matrix', 'perturb', 'perturb1', 'perturb2', 'hybrid'};
            app.OptMethodDropDown.ValueChangedFcn = createCallbackFcn(app, @OptMethodDropDownValueChanged, true);
            app.OptMethodDropDown.BackgroundColor = [0.9412 0.9412 0.9412];
            app.OptMethodDropDown.Layout.Row = 1;
            app.OptMethodDropDown.Layout.Column = 2;
            app.OptMethodDropDown.Value = 'perturb2';

            % Create OptnKnotsEditFieldLabel
            app.OptnKnotsEditFieldLabel = uilabel(app.GridLayoutOpt);
            app.OptnKnotsEditFieldLabel.HorizontalAlignment = 'center';
            app.OptnKnotsEditFieldLabel.Layout.Row = 2;
            app.OptnKnotsEditFieldLabel.Layout.Column = 1;
            app.OptnKnotsEditFieldLabel.Text = 'Opt.nKnots';

            % Create OptnKnotsEditField
            app.OptnKnotsEditField = uieditfield(app.GridLayoutOpt, 'numeric');
            app.OptnKnotsEditField.Limits = [0 1023];
            app.OptnKnotsEditField.RoundFractionalValues = 'on';
            app.OptnKnotsEditField.ValueDisplayFormat = '%.0f';
            app.OptnKnotsEditField.ValueChangedFcn = createCallbackFcn(app, @OptnKnotsEditFieldValueChanged, true);
            app.OptnKnotsEditField.HorizontalAlignment = 'center';
            app.OptnKnotsEditField.Layout.Row = 2;
            app.OptnKnotsEditField.Layout.Column = 2;
            app.OptnKnotsEditField.Value = 23;

            % Create CalculationmethodLabel
            app.CalculationmethodLabel = uilabel(app.GridLayoutOpt);
            app.CalculationmethodLabel.Layout.Row = 1;
            app.CalculationmethodLabel.Layout.Column = 3;
            app.CalculationmethodLabel.Text = 'Calculation method';

            % Create NumberofknotsLabel
            app.NumberofknotsLabel = uilabel(app.GridLayoutOpt);
            app.NumberofknotsLabel.Layout.Row = 2;
            app.NumberofknotsLabel.Layout.Column = 3;
            app.NumberofknotsLabel.Text = 'Number of knots';

            % Create SpinSystemPanel
            app.SpinSystemPanel = uipanel(app.GridLayoutLeftPanel);
            app.SpinSystemPanel.ForegroundColor = [0.502 0.502 0.502];
            app.SpinSystemPanel.Title = 'Spin System';
            app.SpinSystemPanel.BackgroundColor = [1 1 1];
            app.SpinSystemPanel.Layout.Row = 1;
            app.SpinSystemPanel.Layout.Column = 1;
            app.SpinSystemPanel.FontWeight = 'bold';

            % Create GridLayoutSpinSystem
            app.GridLayoutSpinSystem = uigridlayout(app.SpinSystemPanel);
            app.GridLayoutSpinSystem.ColumnWidth = {'2.5x', '3x', '3x'};
            app.GridLayoutSpinSystem.RowHeight = {'1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x'};
            app.GridLayoutSpinSystem.RowSpacing = 5;
            app.GridLayoutSpinSystem.Padding = [2 2 2 2];

            % Create ElectronspinLabel
            app.ElectronspinLabel = uilabel(app.GridLayoutSpinSystem);
            app.ElectronspinLabel.Layout.Row = 1;
            app.ElectronspinLabel.Layout.Column = 3;
            app.ElectronspinLabel.Text = 'Electron spin';

            % Create ElectronspingtensorLabel
            app.ElectronspingtensorLabel = uilabel(app.GridLayoutSpinSystem);
            app.ElectronspingtensorLabel.Layout.Row = 2;
            app.ElectronspingtensorLabel.Layout.Column = 3;
            app.ElectronspingtensorLabel.Text = 'Electron spin g-tensor';

            % Create SysgEditFieldLabel
            app.SysgEditFieldLabel = uilabel(app.GridLayoutSpinSystem);
            app.SysgEditFieldLabel.HorizontalAlignment = 'center';
            app.SysgEditFieldLabel.Layout.Row = 2;
            app.SysgEditFieldLabel.Layout.Column = 1;
            app.SysgEditFieldLabel.Text = 'Sys.g';

            % Create SysgEditField
            app.SysgEditField = uieditfield(app.GridLayoutSpinSystem, 'text');
            app.SysgEditField.ValueChangedFcn = createCallbackFcn(app, @SysgEditFieldValueChanged, true);
            app.SysgEditField.HorizontalAlignment = 'center';
            app.SysgEditField.Layout.Row = 2;
            app.SysgEditField.Layout.Column = 2;
            app.SysgEditField.Value = '2';

            % Create SysNucsEditFieldLabel
            app.SysNucsEditFieldLabel = uilabel(app.GridLayoutSpinSystem);
            app.SysNucsEditFieldLabel.HorizontalAlignment = 'center';
            app.SysNucsEditFieldLabel.Layout.Row = 3;
            app.SysNucsEditFieldLabel.Layout.Column = 1;
            app.SysNucsEditFieldLabel.Text = 'Sys.Nucs';

            % Create SysNucsEditField
            app.SysNucsEditField = uieditfield(app.GridLayoutSpinSystem, 'text');
            app.SysNucsEditField.ValueChangedFcn = createCallbackFcn(app, @SysNucsEditFieldValueChanged, true);
            app.SysNucsEditField.HorizontalAlignment = 'center';
            app.SysNucsEditField.Layout.Row = 3;
            app.SysNucsEditField.Layout.Column = 2;

            % Create SysnEditFieldLabel
            app.SysnEditFieldLabel = uilabel(app.GridLayoutSpinSystem);
            app.SysnEditFieldLabel.HorizontalAlignment = 'center';
            app.SysnEditFieldLabel.Layout.Row = 4;
            app.SysnEditFieldLabel.Layout.Column = 1;
            app.SysnEditFieldLabel.Text = 'Sys.n';

            % Create SysnEditField
            app.SysnEditField = uieditfield(app.GridLayoutSpinSystem, 'text');
            app.SysnEditField.ValueChangedFcn = createCallbackFcn(app, @SysnEditFieldValueChanged, true);
            app.SysnEditField.HorizontalAlignment = 'center';
            app.SysnEditField.Layout.Row = 4;
            app.SysnEditField.Layout.Column = 2;

            % Create NucleiLabel
            app.NucleiLabel = uilabel(app.GridLayoutSpinSystem);
            app.NucleiLabel.Layout.Row = 3;
            app.NucleiLabel.Layout.Column = 3;
            app.NucleiLabel.Text = 'Nuclei';

            % Create ofequivalentnucleiLabel
            app.ofequivalentnucleiLabel = uilabel(app.GridLayoutSpinSystem);
            app.ofequivalentnucleiLabel.Layout.Row = 4;
            app.ofequivalentnucleiLabel.Layout.Column = 3;
            app.ofequivalentnucleiLabel.Text = '# of equivalent nuclei';

            % Create NuclearspinLabel
            app.NuclearspinLabel = uilabel(app.GridLayoutSpinSystem);
            app.NuclearspinLabel.Layout.Row = 5;
            app.NuclearspinLabel.Layout.Column = 3;
            app.NuclearspinLabel.Text = 'Nuclear spin';

            % Create NucleargvalueLabel
            app.NucleargvalueLabel = uilabel(app.GridLayoutSpinSystem);
            app.NucleargvalueLabel.Layout.Row = 6;
            app.NucleargvalueLabel.Layout.Column = 3;
            app.NucleargvalueLabel.Text = 'Nuclear g-value';

            % Create SysAEditFieldLabel
            app.SysAEditFieldLabel = uilabel(app.GridLayoutSpinSystem);
            app.SysAEditFieldLabel.HorizontalAlignment = 'center';
            app.SysAEditFieldLabel.Layout.Row = 7;
            app.SysAEditFieldLabel.Layout.Column = 1;
            app.SysAEditFieldLabel.Text = 'Sys.A';

            % Create SysAEditField
            app.SysAEditField = uieditfield(app.GridLayoutSpinSystem, 'text');
            app.SysAEditField.ValueChangedFcn = createCallbackFcn(app, @SysAEditFieldValueChanged, true);
            app.SysAEditField.HorizontalAlignment = 'center';
            app.SysAEditField.Layout.Row = 7;
            app.SysAEditField.Layout.Column = 2;

            % Create SysQEditFieldLabel
            app.SysQEditFieldLabel = uilabel(app.GridLayoutSpinSystem);
            app.SysQEditFieldLabel.HorizontalAlignment = 'center';
            app.SysQEditFieldLabel.Layout.Row = 8;
            app.SysQEditFieldLabel.Layout.Column = 1;
            app.SysQEditFieldLabel.Text = 'Sys.Q';

            % Create SysQEditField
            app.SysQEditField = uieditfield(app.GridLayoutSpinSystem, 'text');
            app.SysQEditField.ValueChangedFcn = createCallbackFcn(app, @SysQEditFieldValueChanged, true);
            app.SysQEditField.HorizontalAlignment = 'center';
            app.SysQEditField.Layout.Row = 8;
            app.SysQEditField.Layout.Column = 2;

            % Create SysDEditFieldLabel
            app.SysDEditFieldLabel = uilabel(app.GridLayoutSpinSystem);
            app.SysDEditFieldLabel.HorizontalAlignment = 'center';
            app.SysDEditFieldLabel.Layout.Row = 9;
            app.SysDEditFieldLabel.Layout.Column = 1;
            app.SysDEditFieldLabel.Text = 'Sys.D';

            % Create SysDEditField
            app.SysDEditField = uieditfield(app.GridLayoutSpinSystem, 'text');
            app.SysDEditField.ValueChangedFcn = createCallbackFcn(app, @SysDEditFieldValueChanged, true);
            app.SysDEditField.HorizontalAlignment = 'center';
            app.SysDEditField.Layout.Row = 9;
            app.SysDEditField.Layout.Column = 2;

            % Create HyperfinetensorMHzLabel
            app.HyperfinetensorMHzLabel = uilabel(app.GridLayoutSpinSystem);
            app.HyperfinetensorMHzLabel.Layout.Row = 7;
            app.HyperfinetensorMHzLabel.Layout.Column = 3;
            app.HyperfinetensorMHzLabel.Text = 'Hyperfine tensor / MHz';

            % Create QuadrupoletensorMHzLabel
            app.QuadrupoletensorMHzLabel = uilabel(app.GridLayoutSpinSystem);
            app.QuadrupoletensorMHzLabel.Layout.Row = 8;
            app.QuadrupoletensorMHzLabel.Layout.Column = 3;
            app.QuadrupoletensorMHzLabel.Text = 'Quadrupole tensor / MHz';

            % Create ZFScm1Label
            app.ZFScm1Label = uilabel(app.GridLayoutSpinSystem);
            app.ZFScm1Label.Layout.Row = 9;
            app.ZFScm1Label.Layout.Column = 3;
            app.ZFScm1Label.Text = 'ZFS / cm-1';

            % Create LinewidthmTLabel
            app.LinewidthmTLabel = uilabel(app.GridLayoutSpinSystem);
            app.LinewidthmTLabel.Layout.Row = 10;
            app.LinewidthmTLabel.Layout.Column = 3;
            app.LinewidthmTLabel.Text = 'Line width / mT';

            % Create GridLayoutNucSpin
            app.GridLayoutNucSpin = uigridlayout(app.GridLayoutSpinSystem);
            app.GridLayoutNucSpin.RowHeight = {'1x'};
            app.GridLayoutNucSpin.Padding = [0 0 0 0];
            app.GridLayoutNucSpin.Layout.Row = 5;
            app.GridLayoutNucSpin.Layout.Column = 2;

            % Create ILabel
            app.ILabel = uilabel(app.GridLayoutNucSpin);
            app.ILabel.HorizontalAlignment = 'right';
            app.ILabel.Layout.Row = 1;
            app.ILabel.Layout.Column = 1;
            app.ILabel.Text = 'I = ';

            % Create NucSpinVal
            app.NucSpinVal = uilabel(app.GridLayoutNucSpin);
            app.NucSpinVal.HorizontalAlignment = 'center';
            app.NucSpinVal.Layout.Row = 1;
            app.NucSpinVal.Layout.Column = 2;
            app.NucSpinVal.Text = '';

            % Create GridLayoutGn
            app.GridLayoutGn = uigridlayout(app.GridLayoutSpinSystem);
            app.GridLayoutGn.RowHeight = {'1x'};
            app.GridLayoutGn.Padding = [0 0 0 0];
            app.GridLayoutGn.Layout.Row = 6;
            app.GridLayoutGn.Layout.Column = 2;

            % Create gnLabel
            app.gnLabel = uilabel(app.GridLayoutGn);
            app.gnLabel.HorizontalAlignment = 'right';
            app.gnLabel.Layout.Row = 1;
            app.gnLabel.Layout.Column = 1;
            app.gnLabel.Text = 'gn = ';

            % Create NucGnVal
            app.NucGnVal = uilabel(app.GridLayoutGn);
            app.NucGnVal.HorizontalAlignment = 'center';
            app.NucGnVal.Layout.Row = 1;
            app.NucGnVal.Layout.Column = 2;
            app.NucGnVal.Text = '';

            % Create SysSEditFieldLabel
            app.SysSEditFieldLabel = uilabel(app.GridLayoutSpinSystem);
            app.SysSEditFieldLabel.HorizontalAlignment = 'center';
            app.SysSEditFieldLabel.Layout.Row = 1;
            app.SysSEditFieldLabel.Layout.Column = 1;
            app.SysSEditFieldLabel.Text = 'Sys.S';

            % Create SysSEditField
            app.SysSEditField = uieditfield(app.GridLayoutSpinSystem, 'text');
            app.SysSEditField.ValueChangedFcn = createCallbackFcn(app, @SysSEditFieldValueChanged, true);
            app.SysSEditField.HorizontalAlignment = 'center';
            app.SysSEditField.Layout.Row = 1;
            app.SysSEditField.Layout.Column = 2;
            app.SysSEditField.Value = '1/2';

            % Create SyslwEditFieldLabel
            app.SyslwEditFieldLabel = uilabel(app.GridLayoutSpinSystem);
            app.SyslwEditFieldLabel.HorizontalAlignment = 'center';
            app.SyslwEditFieldLabel.Layout.Row = 10;
            app.SyslwEditFieldLabel.Layout.Column = 1;
            app.SyslwEditFieldLabel.Text = 'Sys.lw';

            % Create SyslwEditField
            app.SyslwEditField = uieditfield(app.GridLayoutSpinSystem, 'text');
            app.SyslwEditField.ValueChangedFcn = createCallbackFcn(app, @SyslwEditFieldValueChanged, true);
            app.SyslwEditField.HorizontalAlignment = 'center';
            app.SyslwEditField.Layout.Row = 10;
            app.SyslwEditField.Layout.Column = 2;
            app.SyslwEditField.Value = '0.5';

            % Create ExperimentalParametersPanel
            app.ExperimentalParametersPanel = uipanel(app.GridLayoutLeftPanel);
            app.ExperimentalParametersPanel.ForegroundColor = [0.502 0.502 0.502];
            app.ExperimentalParametersPanel.Title = 'Experimental Parameters';
            app.ExperimentalParametersPanel.BackgroundColor = [1 1 1];
            app.ExperimentalParametersPanel.Layout.Row = 2;
            app.ExperimentalParametersPanel.Layout.Column = 1;
            app.ExperimentalParametersPanel.FontWeight = 'bold';

            % Create GridLayoutExp
            app.GridLayoutExp = uigridlayout(app.ExperimentalParametersPanel);
            app.GridLayoutExp.ColumnWidth = {'2.5x', '3x', '3x'};
            app.GridLayoutExp.RowHeight = {'1x', '1x', '1x', '1x', '1x'};
            app.GridLayoutExp.RowSpacing = 5;
            app.GridLayoutExp.Padding = [2 2 2 2];

            % Create ExpmwFreqEditFieldLabel
            app.ExpmwFreqEditFieldLabel = uilabel(app.GridLayoutExp);
            app.ExpmwFreqEditFieldLabel.HorizontalAlignment = 'center';
            app.ExpmwFreqEditFieldLabel.Layout.Row = 1;
            app.ExpmwFreqEditFieldLabel.Layout.Column = 1;
            app.ExpmwFreqEditFieldLabel.Text = 'Exp.mwFreq';

            % Create ExpmwFreqEditField
            app.ExpmwFreqEditField = uieditfield(app.GridLayoutExp, 'numeric');
            app.ExpmwFreqEditField.Limits = [0 2000];
            app.ExpmwFreqEditField.ValueDisplayFormat = '%11.6g';
            app.ExpmwFreqEditField.ValueChangedFcn = createCallbackFcn(app, @ExpmwFreqEditFieldValueChanged, true);
            app.ExpmwFreqEditField.HorizontalAlignment = 'center';
            app.ExpmwFreqEditField.Layout.Row = 1;
            app.ExpmwFreqEditField.Layout.Column = 2;
            app.ExpmwFreqEditField.Value = 9.4;

            % Create ExpRangeEditFieldLabel
            app.ExpRangeEditFieldLabel = uilabel(app.GridLayoutExp);
            app.ExpRangeEditFieldLabel.HorizontalAlignment = 'center';
            app.ExpRangeEditFieldLabel.Layout.Row = 2;
            app.ExpRangeEditFieldLabel.Layout.Column = 1;
            app.ExpRangeEditFieldLabel.Text = 'Exp.Range';

            % Create ExpRangeEditField
            app.ExpRangeEditField = uieditfield(app.GridLayoutExp, 'text');
            app.ExpRangeEditField.ValueChangedFcn = createCallbackFcn(app, @ExpRangeEditFieldValueChanged, true);
            app.ExpRangeEditField.HorizontalAlignment = 'center';
            app.ExpRangeEditField.Layout.Row = 2;
            app.ExpRangeEditField.Layout.Column = 2;
            app.ExpRangeEditField.Value = '[320 340]';

            % Create mwfrequencyGHzLabel
            app.mwfrequencyGHzLabel = uilabel(app.GridLayoutExp);
            app.mwfrequencyGHzLabel.Layout.Row = 1;
            app.mwfrequencyGHzLabel.Layout.Column = 3;
            app.mwfrequencyGHzLabel.Text = 'm.w. frequency / GHz';

            % Create FieldrangemTLabel
            app.FieldrangemTLabel = uilabel(app.GridLayoutExp);
            app.FieldrangemTLabel.Layout.Row = 2;
            app.FieldrangemTLabel.Layout.Column = 3;
            app.FieldrangemTLabel.Text = 'Field range / mT';

            % Create ExpnPointsEditFieldLabel
            app.ExpnPointsEditFieldLabel = uilabel(app.GridLayoutExp);
            app.ExpnPointsEditFieldLabel.HorizontalAlignment = 'center';
            app.ExpnPointsEditFieldLabel.Layout.Row = 3;
            app.ExpnPointsEditFieldLabel.Layout.Column = 1;
            app.ExpnPointsEditFieldLabel.Text = 'Exp.nPoints';

            % Create ExpnPointsEditField
            app.ExpnPointsEditField = uieditfield(app.GridLayoutExp, 'numeric');
            app.ExpnPointsEditField.Limits = [0 10000];
            app.ExpnPointsEditField.RoundFractionalValues = 'on';
            app.ExpnPointsEditField.ValueDisplayFormat = '%.0f';
            app.ExpnPointsEditField.ValueChangedFcn = createCallbackFcn(app, @ExpnPointsEditFieldValueChanged, true);
            app.ExpnPointsEditField.HorizontalAlignment = 'center';
            app.ExpnPointsEditField.Layout.Row = 3;
            app.ExpnPointsEditField.Layout.Column = 2;
            app.ExpnPointsEditField.Value = 501;

            % Create ExpTemperatureEditFieldLabel
            app.ExpTemperatureEditFieldLabel = uilabel(app.GridLayoutExp);
            app.ExpTemperatureEditFieldLabel.HorizontalAlignment = 'center';
            app.ExpTemperatureEditFieldLabel.Layout.Row = 4;
            app.ExpTemperatureEditFieldLabel.Layout.Column = 1;
            app.ExpTemperatureEditFieldLabel.Text = 'Exp.Temperature';

            % Create ExpTemperatureEditField
            app.ExpTemperatureEditField = uieditfield(app.GridLayoutExp, 'text');
            app.ExpTemperatureEditField.ValueChangedFcn = createCallbackFcn(app, @ExpTemperatureEditFieldValueChanged, true);
            app.ExpTemperatureEditField.HorizontalAlignment = 'center';
            app.ExpTemperatureEditField.Layout.Row = 4;
            app.ExpTemperatureEditField.Layout.Column = 2;

            % Create ExpCrystalOrientationEditFieldLabel
            app.ExpCrystalOrientationEditFieldLabel = uilabel(app.GridLayoutExp);
            app.ExpCrystalOrientationEditFieldLabel.HorizontalAlignment = 'center';
            app.ExpCrystalOrientationEditFieldLabel.Layout.Row = 5;
            app.ExpCrystalOrientationEditFieldLabel.Layout.Column = 1;
            app.ExpCrystalOrientationEditFieldLabel.Text = 'Exp.CrystalOrientation';

            % Create ExpCrystalOrientationEditField
            app.ExpCrystalOrientationEditField = uieditfield(app.GridLayoutExp, 'text');
            app.ExpCrystalOrientationEditField.ValueChangedFcn = createCallbackFcn(app, @ExpCrystalOrientationEditFieldValueChanged, true);
            app.ExpCrystalOrientationEditField.HorizontalAlignment = 'center';
            app.ExpCrystalOrientationEditField.Enable = 'off';
            app.ExpCrystalOrientationEditField.Layout.Row = 5;
            app.ExpCrystalOrientationEditField.Layout.Column = 2;
            app.ExpCrystalOrientationEditField.Value = '[0 0 0]';

            % Create NumberofpointsLabel
            app.NumberofpointsLabel = uilabel(app.GridLayoutExp);
            app.NumberofpointsLabel.Layout.Row = 3;
            app.NumberofpointsLabel.Layout.Column = 3;
            app.NumberofpointsLabel.Text = 'Number of points';

            % Create TemperatureKLabel
            app.TemperatureKLabel = uilabel(app.GridLayoutExp);
            app.TemperatureKLabel.Layout.Row = 4;
            app.TemperatureKLabel.Layout.Column = 3;
            app.TemperatureKLabel.Text = 'Temperature / K';

            % Create EuleranglesdegLabel
            app.EuleranglesdegLabel = uilabel(app.GridLayoutExp);
            app.EuleranglesdegLabel.Layout.Row = 5;
            app.EuleranglesdegLabel.Layout.Column = 3;
            app.EuleranglesdegLabel.Text = 'Euler angles / deg';

            % Create MessageField
            app.MessageField = uieditfield(app.GridLayoutLeftPanel, 'text');
            app.MessageField.Editable = 'off';
            app.MessageField.FontSize = 14;
            app.MessageField.FontColor = [0.149 0.149 0.149];
            app.MessageField.Layout.Row = 6;
            app.MessageField.Layout.Column = 1;
            app.MessageField.Value = '>>';

            % Create GridLayout7
            app.GridLayout7 = uigridlayout(app.GridLayoutLeftPanel);
            app.GridLayout7.RowHeight = {'1x'};
            app.GridLayout7.Padding = [0 0 0 0];
            app.GridLayout7.Layout.Row = 4;
            app.GridLayout7.Layout.Column = 1;

            % Create FrequencyBandPanel
            app.FrequencyBandPanel = uipanel(app.GridLayout7);
            app.FrequencyBandPanel.ForegroundColor = [0.502 0.502 0.502];
            app.FrequencyBandPanel.Title = 'Frequency Band';
            app.FrequencyBandPanel.BackgroundColor = [1 1 1];
            app.FrequencyBandPanel.Layout.Row = 1;
            app.FrequencyBandPanel.Layout.Column = 1;
            app.FrequencyBandPanel.FontWeight = 'bold';

            % Create GridLayout8
            app.GridLayout8 = uigridlayout(app.FrequencyBandPanel);
            app.GridLayout8.ColumnWidth = {'1x'};
            app.GridLayout8.RowHeight = {'1x'};
            app.GridLayout8.Padding = [10 2 10 2];

            % Create FrequencyBandDropDown
            app.FrequencyBandDropDown = uidropdown(app.GridLayout8);
            app.FrequencyBandDropDown.Items = {'S-band', 'X-band', 'Q-band', 'W-band', '240 GHz'};
            app.FrequencyBandDropDown.ValueChangedFcn = createCallbackFcn(app, @FrequencyBandDropDownValueChanged, true);
            app.FrequencyBandDropDown.Layout.Row = 1;
            app.FrequencyBandDropDown.Layout.Column = 1;
            app.FrequencyBandDropDown.Value = 'X-band';

            % Create ExamplesPanel
            app.ExamplesPanel = uipanel(app.GridLayout7);
            app.ExamplesPanel.ForegroundColor = [0.502 0.502 0.502];
            app.ExamplesPanel.Title = 'Examples';
            app.ExamplesPanel.BackgroundColor = [1 1 1];
            app.ExamplesPanel.Layout.Row = 1;
            app.ExamplesPanel.Layout.Column = 2;
            app.ExamplesPanel.FontWeight = 'bold';

            % Create GridLayout9
            app.GridLayout9 = uigridlayout(app.ExamplesPanel);
            app.GridLayout9.ColumnWidth = {'1x'};
            app.GridLayout9.RowHeight = {'1x'};
            app.GridLayout9.Padding = [10 2 10 2];

            % Create ExamplesDropDown
            app.ExamplesDropDown = uidropdown(app.GridLayout9);
            app.ExamplesDropDown.Items = {'Free electron', '1 Proton', '2 Protons', 'Nitroxide radical', 'Mn(III) ion', 'Fe(III) ion', 'Spin Triplet'};
            app.ExamplesDropDown.ValueChangedFcn = createCallbackFcn(app, @ExamplesDropDownValueChanged, true);
            app.ExamplesDropDown.Layout.Row = 1;
            app.ExamplesDropDown.Layout.Column = 1;
            app.ExamplesDropDown.Value = 'Free electron';

            % Create ControlPanel
            app.ControlPanel = uipanel(app.GridLayoutLeftPanel);
            app.ControlPanel.ForegroundColor = [0.502 0.502 0.502];
            app.ControlPanel.Title = 'Control Panel';
            app.ControlPanel.BackgroundColor = [1 1 1];
            app.ControlPanel.Layout.Row = 5;
            app.ControlPanel.Layout.Column = 1;
            app.ControlPanel.FontWeight = 'bold';

            % Create GridLayout10
            app.GridLayout10 = uigridlayout(app.ControlPanel);
            app.GridLayout10.Padding = [10 2 10 2];

            % Create GridLayout11
            app.GridLayout11 = uigridlayout(app.GridLayout10);
            app.GridLayout11.ColumnWidth = {'2.5x', '1.5x'};
            app.GridLayout11.RowHeight = {'1x'};
            app.GridLayout11.Padding = [0 0 0 0];
            app.GridLayout11.Layout.Row = 1;
            app.GridLayout11.Layout.Column = 1;

            % Create SingleorientationSwitchLabel
            app.SingleorientationSwitchLabel = uilabel(app.GridLayout11);
            app.SingleorientationSwitchLabel.Layout.Row = 1;
            app.SingleorientationSwitchLabel.Layout.Column = 1;
            app.SingleorientationSwitchLabel.Text = 'Single orientation';

            % Create SingleorientationSwitch
            app.SingleorientationSwitch = uiswitch(app.GridLayout11, 'slider');
            app.SingleorientationSwitch.ValueChangedFcn = createCallbackFcn(app, @SingleorientationSwitchValueChanged, true);
            app.SingleorientationSwitch.Layout.Row = 1;
            app.SingleorientationSwitch.Layout.Column = 2;

            % Create GridLayout11_2
            app.GridLayout11_2 = uigridlayout(app.GridLayout10);
            app.GridLayout11_2.ColumnWidth = {'2.5x', '1.5x'};
            app.GridLayout11_2.RowHeight = {'1x'};
            app.GridLayout11_2.Padding = [0 0 0 0];
            app.GridLayout11_2.Layout.Row = 2;
            app.GridLayout11_2.Layout.Column = 1;

            % Create OverlayresultsSwitchLabel
            app.OverlayresultsSwitchLabel = uilabel(app.GridLayout11_2);
            app.OverlayresultsSwitchLabel.Layout.Row = 1;
            app.OverlayresultsSwitchLabel.Layout.Column = 1;
            app.OverlayresultsSwitchLabel.Text = 'Overlay results';

            % Create OverlayresultsSwitch
            app.OverlayresultsSwitch = uiswitch(app.GridLayout11_2, 'slider');
            app.OverlayresultsSwitch.ValueChangedFcn = createCallbackFcn(app, @OverlayresultsSwitchValueChanged, true);
            app.OverlayresultsSwitch.Layout.Row = 1;
            app.OverlayresultsSwitch.Layout.Column = 2;

            % Create GridLayout12
            app.GridLayout12 = uigridlayout(app.GridLayout10);
            app.GridLayout12.ColumnWidth = {'1x'};
            app.GridLayout12.RowHeight = {'1x'};
            app.GridLayout12.Padding = [10 0 0 0];
            app.GridLayout12.Layout.Row = 1;
            app.GridLayout12.Layout.Column = 2;

            % Create StartButton
            app.StartButton = uibutton(app.GridLayout12, 'state');
            app.StartButton.ValueChangedFcn = createCallbackFcn(app, @StartButtonValueChanged, true);
            app.StartButton.Text = 'Start';
            app.StartButton.FontWeight = 'bold';
            app.StartButton.Layout.Row = 1;
            app.StartButton.Layout.Column = 1;

            % Create GridLayout13
            app.GridLayout13 = uigridlayout(app.GridLayout10);
            app.GridLayout13.ColumnWidth = {'8x', '1x', '6x'};
            app.GridLayout13.RowHeight = {'1x'};
            app.GridLayout13.Padding = [10 0 0 0];
            app.GridLayout13.Layout.Row = 2;
            app.GridLayout13.Layout.Column = 2;

            % Create DataisloadedLampLabel
            app.DataisloadedLampLabel = uilabel(app.GridLayout13);
            app.DataisloadedLampLabel.HorizontalAlignment = 'right';
            app.DataisloadedLampLabel.Layout.Row = 1;
            app.DataisloadedLampLabel.Layout.Column = 1;
            app.DataisloadedLampLabel.Text = 'Data is loaded';

            % Create DataisloadedLamp
            app.DataisloadedLamp = uilamp(app.GridLayout13);
            app.DataisloadedLamp.Layout.Row = 1;
            app.DataisloadedLamp.Layout.Column = 2;
            app.DataisloadedLamp.Color = [0.502 0.502 0.502];

            % Create ClearButton
            app.ClearButton = uibutton(app.GridLayout13, 'push');
            app.ClearButton.ButtonPushedFcn = createCallbackFcn(app, @ClearButtonPushed, true);
            app.ClearButton.Layout.Row = 1;
            app.ClearButton.Layout.Column = 3;
            app.ClearButton.Text = 'Clear';

            % Create EnergyleveldiagramPanel
            app.EnergyleveldiagramPanel = uipanel(app.GridLayoutMain);
            app.EnergyleveldiagramPanel.ForegroundColor = [0.502 0.502 0.502];
            app.EnergyleveldiagramPanel.Title = 'Energy level diagram';
            app.EnergyleveldiagramPanel.BackgroundColor = [1 1 1];
            app.EnergyleveldiagramPanel.Layout.Row = 1;
            app.EnergyleveldiagramPanel.Layout.Column = 2;
            app.EnergyleveldiagramPanel.FontWeight = 'bold';

            % Create GridLayoutLevels
            app.GridLayoutLevels = uigridlayout(app.EnergyleveldiagramPanel);
            app.GridLayoutLevels.ColumnWidth = {'1x'};
            app.GridLayoutLevels.RowHeight = {'1x'};

            % Create UIAxesLevels
            app.UIAxesLevels = uiaxes(app.GridLayoutLevels);
            title(app.UIAxesLevels, 'Selected orientation')
            xlabel(app.UIAxesLevels, 'Magnetic field / mT')
            ylabel(app.UIAxesLevels, 'Energy / cm^{-1}')
            app.UIAxesLevels.XLim = [320 340];
            app.UIAxesLevels.Box = 'on';
            app.UIAxesLevels.Layout.Row = 1;
            app.UIAxesLevels.Layout.Column = 1;

            % Create GridLayoutRightCol
            app.GridLayoutRightCol = uigridlayout(app.GridLayoutMain);
            app.GridLayoutRightCol.ColumnWidth = {'1x'};
            app.GridLayoutRightCol.Padding = [0 0 0 0];
            app.GridLayoutRightCol.Layout.Row = 1;
            app.GridLayoutRightCol.Layout.Column = 3;

            % Create AbsorptionspectrumPanel
            app.AbsorptionspectrumPanel = uipanel(app.GridLayoutRightCol);
            app.AbsorptionspectrumPanel.ForegroundColor = [0.502 0.502 0.502];
            app.AbsorptionspectrumPanel.Title = 'Absorption spectrum';
            app.AbsorptionspectrumPanel.BackgroundColor = [1 1 1];
            app.AbsorptionspectrumPanel.Layout.Row = 1;
            app.AbsorptionspectrumPanel.Layout.Column = 1;
            app.AbsorptionspectrumPanel.FontWeight = 'bold';

            % Create GridLayoutAbsorption
            app.GridLayoutAbsorption = uigridlayout(app.AbsorptionspectrumPanel);
            app.GridLayoutAbsorption.ColumnWidth = {'1x'};
            app.GridLayoutAbsorption.RowHeight = {'1x'};

            % Create UIAxesAbsorption
            app.UIAxesAbsorption = uiaxes(app.GridLayoutAbsorption);
            title(app.UIAxesAbsorption, 'Powder spectrum')
            xlabel(app.UIAxesAbsorption, 'Magnetic field / mT')
            app.UIAxesAbsorption.PlotBoxAspectRatio = [1.09126980177839 1 1];
            app.UIAxesAbsorption.XLim = [320 340];
            app.UIAxesAbsorption.YLim = [-0.2 1.2];
            app.UIAxesAbsorption.Box = 'on';
            app.UIAxesAbsorption.Layout.Row = 1;
            app.UIAxesAbsorption.Layout.Column = 1;

            % Create CWspectrumPanel
            app.CWspectrumPanel = uipanel(app.GridLayoutRightCol);
            app.CWspectrumPanel.ForegroundColor = [0.502 0.502 0.502];
            app.CWspectrumPanel.Title = 'CW spectrum';
            app.CWspectrumPanel.BackgroundColor = [1 1 1];
            app.CWspectrumPanel.Layout.Row = 2;
            app.CWspectrumPanel.Layout.Column = 1;
            app.CWspectrumPanel.FontWeight = 'bold';
            app.CWspectrumPanel.FontSize = 14;

            % Create GridLayoutCW
            app.GridLayoutCW = uigridlayout(app.CWspectrumPanel);
            app.GridLayoutCW.ColumnWidth = {'1x'};
            app.GridLayoutCW.RowHeight = {'1x'};

            % Create UIAxesCW
            app.UIAxesCW = uiaxes(app.GridLayoutCW);
            title(app.UIAxesCW, 'Powder spectrum')
            xlabel(app.UIAxesCW, 'Magnetic field / mT')
            app.UIAxesCW.PlotBoxAspectRatio = [1.09999995628993 1 1];
            app.UIAxesCW.XLim = [320 340];
            app.UIAxesCW.YLim = [-1.2 1.2];
            app.UIAxesCW.Box = 'on';
            app.UIAxesCW.Layout.Row = 1;
            app.UIAxesCW.Layout.Column = 1;

            % Show the figure after all components are created
            app.EPRSimulatorUIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = eprsimulator

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.EPRSimulatorUIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.EPRSimulatorUIFigure)
        end
    end
end

