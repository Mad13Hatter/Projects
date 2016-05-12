function [t, h] = call_dstate()
 
    time = (0:0.01:400)
    tspan = [time]; % time interval from 0 to 15 seconds
    h0 = 0.15; % initial condition
    
    % dstate evaluates the R.H.S of the ODE in the format (dh/dt)= R.H.S
    [t, h] = ode45(@dstate,tspan,h0); % @dstate is referring to the function dstate (is a function handle)
    [t2, h2] = ode45(@dstate2,tspan,h0); 
    [t3, h3] = ode45(@dstate3,tspan,h0); 
    disp([t, h])
    disp([t2,h2])
    disp([t3,h3])
    close all;
   
    figure(10)
    hold on;
    grid on;
    plot(t,h,'LineWidth',3);
    plot(t,h2,'LineWidth',3);
    xlabel('Time [t] (s)');
    ylabel('Height [h(t)] (s)');
    title('Original h(t) & h(t) after neglecting flow acceleration ');
    
    figure(11)
    hold on;
    grid on;
    plot(t,h,'LineWidth',3);
    plot(t,h3,'-g','LineWidth',3);
    xlabel('Time [t] (s)');
    ylabel('Height [h(t)] (s)');
    title('Original h(t) & h(t) after neglecting entry friction ');
    
        function dhdt = dstate(t,h)
            
            % Parameters
            lpip = (0.01:0.024:0.25)'
            g = 9.8; % [m/s^2]
            h0 = 0.15; % initial height in [m]
            H3 = 0.044; % height from the bottom of the bottle, where the tube in located in [m]
            K1 = 164;
            u = 1.003 * 10 ^ (-3); % viscosity in [N-s/m^2]
            p = 998; % density in [kg/m^3] at room temp
            K2 = 1; % for a re-entrant tube
            e = 0.0000015; % roughness factor for plastic from table 6.1 in [m]
            D = 0.109;
 
            % Variables
            L = 0.04; % varies between 0.01 and 0.25 [m]
            d = 0.004; % varies between 0.004 and 0.01 [m]; does not have to vary, but it is good to show the involvement of this variable with laminar vs turbulent
 
            % ODE
            A = ((D^2)/(p*d^4))*(K1*u*d + 64*u*L);
            B = (D^4/((p^2)*(d^8)))*((K1*u*d + 64*u*L)^2);
            C = 8*(D^4)*((K2 + 1)/(d^4))*g*(H3 - h);
            G = (2*(K2+1)*(D^4/d^4));
            
            dhdt = (A - sqrt(B - C))/G;
        end
    
        function dhdt2 = dstate2(t,h)
            
            % Parameters
            lpip = (0.01:0.024:0.25)'
            g = 9.8; % [m/s^2]
            h0 = 0.15; % initial height in [m]
            H3 = 0.044; % height from the bottom of the bottle, where the tube in located in [m]
            K1 = 164;
            u = 1.003 * 10 ^ (-3); % viscosity in [N-s/m^2]
            p = 998; % density in [kg/m^3] at room temp
            K2 = 1; % for a re-entrant tube
            e = 0.0000015; % roughness factor for plastic from table 6.1 in [m]
            D = 0.109;
 
            % Variables
            L = 0.04; % varies between 0.01 and 0.25 [m]
            d = 0.004; % varies between 0.004 and 0.01 [m]; does not have to vary, but it is good to show the involvement of this variable with laminar vs turbulent
            
            % Ignnoring the Flow acceleration term in Bernoullis equation
            A2 = ((D^2)/(p*d^4))*(K1*u*d + 64*u*L);
            B2 = (D^4/((p^2)*(d^8)))*((K1*u*d + 64*u*L)^2);
            C2 = 8*(D^4)*((K2)/(d^4))*g*(H3 - h);
            G2 = (2*(K2)*(D^4/d^4));
 
            dhdt2 = (A2 - sqrt(B2 - C2))/G2;
        end
        function dhdt3 = dstate3(t,h)
            
            % Parameters
            lpip = (0.01:0.024:0.25)'
            g = 9.8; % [m/s^2]
            h0 = 0.15; % initial height in [m]
            H3 = 0.044; % height from the bottom of the bottle, where the tube in located in [m]
            K1 = 164;
            u = 1.003 * 10 ^ (-3); % viscosity in [N-s/m^2]
            p = 998; % density in [kg/m^3] at room temp
            K2 = 1; % for a re-entrant tube
            e = 0.0000015; % roughness factor for plastic from table 6.1 in [m]
            D = 0.109;
 
            % Variables
            L = 0.04; % varies between 0.01 and 0.25 [m]
            d = 0.004; % varies between 0.004 and 0.01 [m]; does not have to vary, but it is good to show the involvement of this variable with laminar vs turbulent
            
           % Ignnoring the Entry Friction term in Bernoullis equation
            K1 = 0;
            K2 = 0;
            A3 = ((D^2)/(p*d^4))*(K1*u*d + 64*u*L);
            B3 = (D^4/((p^2)*(d^8)))*((K1*u*d + 64*u*L)^2);
            C3 = 8*(D^4)*((K2 + 1)/(d^4))*g*(H3 - h);
            G3 = (2*(K2+1)*(D^4/d^4));
 
            dhdt3 = (A3 - sqrt(B3 - C3))/G3;
        end
end


