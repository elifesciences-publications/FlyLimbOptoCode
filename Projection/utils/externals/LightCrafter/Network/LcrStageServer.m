classdef LcrStageServer < StageServer
    
    properties (Access = private)
        lightCrafter
        background
    end
    
    methods (Static)
        
        function server = createAndStart(monitorNumber, port)
            if nargin < 1
                monitorNumber = 2;
            end
            
            if nargin < 2
                port = 5678;
            end
                        
            server = LcrStageServer(port);
            server.start(Lcr4500.NATIVE_RESOLUTION, true, LcrMonitor(monitorNumber));
        end
        
    end
    
    methods
        
        function obj = LcrStageServer(port)
            if nargin < 1
                port = 5678;
            end
            
            obj = obj@StageServer(port);
        end
        
    end
    
    methods (Access = protected)
        
        function willStart(obj)
            willStart@StageServer(obj);
            
            monitor = obj.canvas.window.monitor;
            monitor.setGamma(1);
            
            obj.lightCrafter = Lcr4500(monitor);
            obj.lightCrafter.connect();
            
            % Set LEDs to enable automatically.
            obj.lightCrafter.setLedEnables(true, true, true, true);
            
            obj.lightCrafter.setMode(LcrMode.PATTERN);
            obj.lightCrafter.setPatternAttributes(Lcr4500.MAX_PATTERN_BIT_DEPTH, 'white', 1);
            
            obj.background = Rectangle();
            obj.background.position = obj.canvas.size/2;
            obj.background.size = obj.canvas.size;
            obj.background.color = 0;
            obj.background.init(obj.canvas);
            
            if monitor.resolution == Lcr4500.NATIVE_RESOLUTION
                % Stretch the projection matrix to account for the LightCrafter diamond pixel screen.
                window = obj.canvas.window;
                obj.canvas.projection.setIdentity();
                obj.canvas.projection.orthographic(0, window.size(1)*2, 0, window.size(2));
                
                obj.background.position(1) = obj.background.position(1)*2;
                obj.background.size(1) = obj.background.size(1)*2;
            end
        end
        
        function didStop(obj)
            didStop@StageServer(obj);
            
            obj.lightCrafter.disconnect();
        end
        
        function onEventReceived(obj, src, data)
            client = data.client;
            value = data.value;
            
            try
                switch value{1}
                    case LcrNetEvents.GET_LCR_PATTERN_ATTRIBUTES
                        obj.onEventGetLcrPatternAttributes(client, value);
                    case LcrNetEvents.SET_LCR_PATTERN_ATTRIBUTES
                        obj.onEventSetLcrPatternAttributes(client, value);
                    case LcrNetEvents.GET_LCR_LED_CURRENTS
                        obj.onEventGetLcrLedCurrents(client, value);
                    case LcrNetEvents.SET_LCR_LED_CURRENTS
                        obj.onEventSetLcrLedCurrents(client, value);
                    case LcrNetEvents.GET_LCR_LED_ENABLES
                        obj.onEventGetLcrLedEnables(client, value);
                    case LcrNetEvents.SET_LCR_LED_ENABLES
                        obj.onEventSetLcrLedEnables(client, value);
                    otherwise
                        onEventReceived@StageServer(obj, src, data);
                end
            catch x
                client.send(NetEvents.ERROR, x);
            end
        end
        
        function onEventGetLcrPatternAttributes(obj, client, value) %#ok<INUSD>
            [bitDepth, color, numPatterns] = obj.lightCrafter.getPatternAttributes();
            client.send(NetEvents.OK, bitDepth, color, numPatterns);
        end
        
        function onEventSetLcrPatternAttributes(obj, client, value)
            bitDepth = value{2};
            color = value{3};
            numPatterns = value{4};
            
            if isempty(numPatterns)
                numPatterns = obj.lightCrafter.maxNumPatternsForBitDepth(bitDepth);
            end
            
            [cBitDepth, cColor, cNumPatterns] = obj.lightCrafter.getPatternAttributes();
            if bitDepth == cBitDepth && strncmpi(color, cColor, length(color)) && numPatterns == cNumPatterns
                client.send(NetEvents.OK);
                return;
            end
            
            if bitDepth ~= cBitDepth || numPatterns ~= cNumPatterns
                obj.sessionData.player = [];
            end
            
            obj.lightCrafter.setPatternAttributes(bitDepth, color, numPatterns);
            client.send(NetEvents.OK);
        end
        
        function onEventGetLcrLedCurrents(obj, client, value) %#ok<INUSD>
            [red, green, blue] = obj.lightCrafter.getLedCurrents();
            client.send(NetEvents.OK, red, green, blue);
        end
        
        function onEventSetLcrLedCurrents(obj, client, value)
            red = value{2};
            green = value{3};
            blue = value{4};
            
            obj.lightCrafter.setLedCurrents(red, green, blue);
            client.send(NetEvents.OK);
        end
        
        function onEventGetLcrLedEnables(obj, client, value) %#ok<INUSD>
            [auto, red, green, blue] = obj.lightCrafter.getLedEnables();
            client.send(NetEvents.OK, auto, red, green, blue);
        end
        
        function onEventSetLcrLedEnables(obj, client, value)
            auto = value{2};
            red = value{3};
            green = value{4};
            blue = value{5};
            
            obj.lightCrafter.setLedEnables(auto, red, green, blue);
            client.send(NetEvents.OK);
        end
        
        function onEventGetCanvasSize(obj, client, value) %#ok<INUSD>
            size = obj.canvas.size;
            if obj.canvas.window.monitor.resolution == Lcr4500.NATIVE_RESOLUTION
                % Stretch for diamond pixel layout.
                size(1) = size(1) * 2;
            end
            
            client.send(NetEvents.OK, size);
        end
        
        function onEventSetCanvasClearColor(obj, client, value)
            color = value{2};
            
            obj.background.color = color;            
            client.send(NetEvents.OK);
        end
        
        function onEventPlay(obj, client, value)
            presentation = value{2};
            prerender = value{3};
            
            % Add the background to the presentation.
            presentation.insertStimulus(1, obj.background);
            
            if prerender
                obj.sessionData.player = PrerenderedPlayer(presentation);
            else
                obj.sessionData.player = RealtimePlayer(presentation);
            end
            
            [bitDepth, ~, nPatterns] = obj.lightCrafter.getPatternAttributes();
            renderer = LcrPatternRenderer(nPatterns, bitDepth);
                        
            obj.canvas.setRenderer(renderer);
            resetRenderer = onCleanup(@()obj.canvas.resetRenderer());
            
            compositor = LcrPatternCompositor();
            compositor.bindPatternRenderer(renderer);
            
            obj.sessionData.player.setCompositor(compositor);
            
            % Unlock client to allow async operations during play.
            client.send(NetEvents.OK);
            
            try
                obj.sessionData.playInfo = obj.sessionData.player.play(obj.canvas);
            catch x
                obj.sessionData.playInfo = x;
            end
        end
        
    end
    
end