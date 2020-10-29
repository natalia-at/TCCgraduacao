function DP_v8()
%------------------------------------------------------------------------------------------------
% UNIT COMMITMENT BY DYNAMIC PROGRAMMING METHOD
%
% Program versions development:
% version v1    - basic DP algorithm with simple quick dispatch (linear, one segment generation curve)
%                 based on the book:
%                 A. J. Wood and B. F. Wollenberg, Power Generation Operation and Control,
%                 1984, John Wiley, New York
% version v2    - minimum up and down times taken into account: this is based on the following
%                 article (with slight modification since cooling time is not taken into account)
%                   C.Li, R.B.Johnson, and A.J.Svoboda: "A new unit commitment method",
%                   IEEE Transactions on Power Systems, Vol. 12, No. 1, 1997, pp. 113–119, 1997.
% version v3    - generators dispatch based on LINPROG
% version v4    - ramp up and down constraints taken into account, dispatching based on LINPROG
% version v4_b  - same as v4 with either simple quick dispatch or linprog based one
% version v5    - enhanced DP method (keep track of several predecessors, not only one)
%                 Version 5 is based on the article:
%                   W.J.Hobbs, et.al., “An enhanced dynamic programming approach for unit commitment,”
%                   IEEE Trans. Power Syst., Vol. 3, No. 3, pp. 1201-1205, August 1988.
% version v6    - deterministic spinning reserve added; also QUADPROG option included
% version v7    - shut down cost taken into account; start-up costs can now be cold/hot or
%                 exponential.
% version V8    - some basic check up on data availability added on. 
%
% Program can use either:
%       - priority list of generators to be commited
%       - complete enumeration (all possible combinations)
% Priority list is based on the no load cost and incremental heat rate data.
%
% This program, uses a forward DP technique, and does not
% take into account time resolution different from 1h.
%
% Program features are controlled by a set of the switches (flags).
%
% THINGS THAT SHOULD BE IMPROVED/CHANGED IN THE CODE:
% - modify the code so it can deal with any time resolution (not only 1h)
% - since there is a relation between coefficients in linear cost methods and coefficients in 
%   quadratic cost method, it is possible to estimate the former one (if not given) 
%   based on the latter one. This can be done off-line (using for example, least square method)
%  and it would offer more versatility to the program (eg. using quick dispatch or priority list
%  even if only quadratic coefficients are given).
% - ... whatever else pops in your mind.
%
% Author: Vladimir Stanojevic, March 2011
%------------------------------------------------------------------------------------------------
%%
tic                                     % initialize timer
clc
warning off
%------------------------------------------------------------------------------------------------
% read input data from a file
DP_input_data_n

GMIN            = gen_data(:, 2);           % generator min power                           [MW]
GMAX            = gen_data(:, 3);           % generator max. power                          [MW]
GINC            = gen_data(:, 4);           % incremental heat rate                         [BTU/kWh]
GNLC            = gen_data(:, 5);           % generator no load cost                        [£/h]
GSC             = gen_data(:, 6);           % generator start up cost (cold), also BETA     [£]
GFC             = gen_data(:, 7);           % generator fuel cost                           [£/cMBTU]
GMINUP          = gen_data(:, 8);           % generator min. up time                        [h]
GMINDOWN        = gen_data(:, 9);           % generator min. down time                      [h]
GSTATINI        = gen_data(:,10);           % generator initial status (time on/off)        [h]
GSH             = gen_data(:,11);           % generator start up cost (hot), also ALPHA     [£]
GCSTIME         = gen_data(:,12);           % generator cold start time                     [h]
GRAMPUP         = gen_data(:,13);           % generator ramp up rate                        [MW/h]
GRAMPDOWN       = gen_data(:,14);           % generator ramp down rate                      [MW/h]
COEF_A          = gen_data(:,15);           % free term in quadratic-cost function          [£]
COEF_B          = gen_data(:,16);           % linear term in quadratic-cost function        [£/MWh]
COEF_C          = gen_data(:,17);           % 2nd order term in quadratic-cost function     [£/(MW^2)h]
GSDC            = gen_data(:,18);           % generator shut down cost                      [£]
TAU             = gen_data(:,19);           % generator shut up cost exp. coef.             [£]
%------------------------------------------------------------------------------------------------
NG              = size(gen_data,1);         % no. of generators
NT              = size(DEMAND,1);           % number of time periods (hours)
%------------------------------------------------------------------------------------------------

%------------------------------------------------------------------------------------------------
% check up the availability of data for certain cases:
if (DISPATCH_METHOD == 2 | DISPATCH_METHOD == 3) & (any(isnan(GNLC)) | any(isnan(GFC)) | any(isnan(GINC)))
    STR = ['To use linear cost model, you must provide data for NO LOAD COSTS,'...
        'FUEL COSTS and INCREMENTAL COSTS.'];   % there are no data for quick dispatch method,
    msgbox(STR,'DATA AVAILABILITY CHECK FAILURE!','warn');       % write a message
    return
elseif DISPATCH_METHOD == 1 & (any(isnan(COEF_A)) | any(isnan(COEF_B)) | any(isnan(COEF_C)))
    STR = ['To use quadratic cost model, you must provide data for the cost coefficients:,'...
        'COEFF_A (£), COEFF_B (£/MWh) and COEFF_C (£/MW^2h).'];   % there are no data for quick dispatch method,
    msgbox(STR,'DATA AVAILABILITY CHECK FAILURE!','warn');       % write a message
    return
end
if MIN_UP_DOWN_TIME_FLAG == 1 & (any(isnan(GMINUP)) | any(isnan(GMINDOWN)))
    STR = ['To use minimum up and down time constraints, you must provide data for GMINUP'...
        ' and GMINDOWN.'];
    msgbox(STR,'DATA AVAILABILITY CHECK FAILURE!','warn');
    return
end
if RAMP_UP_DOWN_FLAG == 1 & (any(isnan(GRAMPUP)) | any(isnan(GRAMPDOWN)))
    STR = ['To use rump constraints, you must provide data for GRAMPUP '...
        'and GRAMPDOWN.'];
    msgbox(STR,'DATA AVAILABILITY CHECK FAILURE!','warn');
    return
end
if COMPLETE_ENUMERATION_FLAG == 0 & (any(isnan(GNLC)) | any(isnan(GFC)) | any(isnan(GINC)))
    STR = ['To use priority list, you must provide data for NO LOAD COSTS,'...
        'FUEL COSTS and INCREMENTAL COSTS.'];
    msgbox(STR,'DATA AVAILABILITY CHECK FAILURE!','warn');
    return
end
if START_UP_COST_METHOD == 2 & (any(isnan(GSH)) | any(isnan(GCSTIME)) )
    STR = ['To use cold/hot start up cost method, you must provide data for GSH '...
        'and GCSTIME.'];
    msgbox(STR,'DATA AVAILABILITY CHECK FAILURE!','warn');
    return
elseif START_UP_COST_METHOD == 3 & (any(isnan(GSH)) | any(isnan(TAU)) )
    STR = ['To use exponential start up cost method, you must provide data for GSH '...
        'and TAU.'];
    msgbox(STR,'DATA AVAILABILITY CHECK FAILURE!','warn');
    return
end
%------------------------------------------------------------------------------------------------

% enabling / disabling some of the constraints
if MIN_UP_DOWN_TIME_FLAG == 0               % minimum up and down time enabled/disabled
    GMINUP(:)   = 1;                        % if disabled, all min up and down times
    GMINDOWN(:) = 1;                        % are set to 1
end
if RAMP_UP_DOWN_FLAG == 0                   % ramping constraints enabled/disabled
    GRAMPUP(:)   = Inf;                     % if disabled, ramp rates are
    GRAMPDOWN(:) = Inf;                     % set to a large number
end
if COMPLETE_ENUMERATION_FLAG == 0           % use either priority list or...
    [GMAXlst,GMINlst,LIST_STATES,GEN_ORDER] = priority_list(GNLC,GFC,GMAX,GMIN,GINC,NG);
else                                        % ...complete enumeration (consisting of all possible combinations)
    [GMAXlst,GMINlst,LIST_STATES,GEN_ORDER] = complete_enumeration(GNLC,GFC,GMAX,GMIN,GINC,NG);
end
% analizes if the spinning reserva will be considered or not
if RESERVE_FLAG == 1
    if exist('RES_UP','var') ~= 1 | isempty(RES_UP) % if reserve-up vector is not defined or if it is empty
        RES_UP = K_RES_UP * DEMAND;                 % create reserve-up vector in proportion to demand
    end
    if exist('RES_DN','var') ~= 1 | isempty(RES_DN) % if reserve down vector is not defined or if it is empty
        RES_DN = K_RES_DN * DEMAND;                 % create reserve down vector in proportion to demand
    end
else
    RES_UP = zeros(size(DEMAND));                   % if reserve not required,
    RES_DN = zeros(size(DEMAND));                   % set it to zero.
end
if START_UP_COST_METHOD == 3                        % if start-up costs are exponential
    ALPHA = GSH;                                    % then define ALPHA
    BETA  = GSC;                                    % and BETA
else
    ALPHA = NaN*ones(NG,1);                         % otherwise, just define the names for variables
    BETA  = NaN*ones(NG,1);                         % since they will be passed to the functions
end
%------------------------------------------------------------------------------------------------

% Determines the initial status (ON/OFF = 1/0) for each generator, based on the input data (GSTATINI).
% (GSTATINI contains the number of hours that a generator was ON/OFF before the 1st time step)
% INI_STATE [NG x 1] - initial states of generators (1-commited, 0-not commited)
% INI_STATE_NUM      - position (column) of vector INI_STATE in the list of states
INI_STATE = (GSTATINI > 0);
[I, INI_STATE_NUM]= ismember(INI_STATE',LIST_STATES','rows');

%------------------------------------------------------------------------------------------------
% main loop of the program - search for optimal commitment by dynamic programming
%------------------------------------------------------------------------------------------------
% Here is the brief explanation of the algorithm:
% State is a unique combination of commited and non-commited generators.
% Commited generators are designated with "1" and non-commited generators are designated with "0".
%
% For each hour, program finds the potentially feasible states. Potentially feasible states are
% the states where demand (and reserve) can be supplied by the commited generators.
% If there are no potentially feasible states, program displays the error message and terminates.
%
% For each potentially feasible state, program takes all feasible states from the previous hour
% and checks if the transition to the current state (in current hour) is possible.
% If it is not possible, the corresponding transition (start-up) cost is set to Inf.
% However, if the transition is possible, calculated are the transition costs. Production for
% the current hour is calculated based on demand taking into account production at previous hour (ramp-up and
% down constraints). Finally, total cost is the sum of the transition cost, production cost, and
% the total cost at the state in previous hour. This procedure is repeated for all the states in
% previous hour. Total costs are then sorted and MN of them are saved (this is enhancement comparing
% to the classical dynamic program where only 1 previous state is saved). If the transition to
% a state in current hour is not possible from any of the states in previous hour, then current state is
% regarded as infeasible and is not considered anymore. If all the states in an hour are infeasible,
% program displays the error message and terminates.
%------------------------------------------------------------------------------------------------
for HOUR = 1:NT
    fprintf('Currently processing hour: %2d \n',HOUR)
    if HOUR == 1
        PREV_STATES_NUM = INI_STATE_NUM;            % Positions (columns) of feasible states in the list of states, for previous hour
        X_PREV  = GSTATINI;                         % number of hours generators are ON (>0) or OFF (<0).
        PRODUCTION_PREV = zeros(size(X_PREV));      % generator outputs
        TR_PREV = PREV_STATES_NUM;                  % transition path matrix
        FC_PREV = 0;                                % cumulative cost vector
    else
        X_PREV = X_CURR;                            % keep the number of gen. working hours for each state (at previous hour)
        PRODUCTION_PREV = PRODUCTION_CURR;          % and the gen. outputs for previous hour states
        TR_PREV = TR;                               % rows of matrix TR define the transition path
        FC_PREV = FC;                               % save the cumulative cost vector from the previous hour
        PREV_STATES_NUM = TR_PREV(1:COUNTER,end);   % states in the previous hour are given in the last column of TR
    end
    % FEASIBLE_STATES_NUM = positions (columns) of potentially feasible states in the list of states, for current hour.
    [FEASIBLE_STATES_NUM,SUCCESS] = find_feasible_states(GMINlst,GMAXlst,DEMAND,HOUR,RES_UP,RES_DN);
    if SUCCESS == 0                                 % if unable to find any feasible state to match demand,
        return                                      % quit the program
    end
    MN = min(length(PREV_STATES_NUM),N_PRED);                      % number of predecessors to be examined
    X_CURR = zeros(NG,MN*length(FEASIBLE_STATES_NUM));             % prepare the number of gen. working hours for each state (at current hour)
    PRODUCTION_CURR = zeros(NG,MN*length(FEASIBLE_STATES_NUM));    % prepare generator outputs for each state at current hour
    TR = zeros(MN*length(FEASIBLE_STATES_NUM),HOUR+1);             % prepare transition path matrix
    FC = zeros(MN*length(FEASIBLE_STATES_NUM),1);                  % prepare cumulative cost vector
    COUNTER = 0;
    % take each feasible (current hour) state and...
    for J = 1: length(FEASIBLE_STATES_NUM)
        GEN_START_SHUT_COST = zeros(NG,1);                         % start up (shut down) costs
        TOTAL_COST = zeros(1,length(PREV_STATES_NUM));             % total cost (production cost + start up cost + total cost of previous hour)
        % X_TEMP - temporarily stores number of gen. working hours for combination of current state and all previous hour states
        X_TEMP = zeros(NG,length(PREV_STATES_NUM));
        % PRODUCTION_TEMP - temporarily stores gen. outputs for combination of current state and all previous hour states
        PRODUCTION_TEMP = zeros(NG,length(PREV_STATES_NUM));
        % take a state (from all feasible states), one by one; let it be CURRENT_STATE
        CURRENT_STATE  = LIST_STATES(:,FEASIBLE_STATES_NUM(J));
        %----------------------------------------------------------------------------------
        % ... compare it with each feasible state at previous hour
        for K = 1: length(PREV_STATES_NUM)
            if HOUR == 1;
                PREVIOUS_STATE = INI_STATE;
            else
                PREVIOUS_STATE = LIST_STATES(:,PREV_STATES_NUM(K));
            end
            % check if the transition from previous state to the current state is possible regarding min up and down times
            [X,SUCCESS] = check_up_down_time(CURRENT_STATE,PREVIOUS_STATE,X_PREV(:,K),GMINUP,GMINDOWN,NG);
            if SUCCESS==0                                                   % if it is not possible,
                GEN_START_SHUT_COST(:,K) = Inf;                             % mark the transition cost and
                PROD_COST = ones(NG,1)*Inf;                                 % production cost as extremely high
                GEN_PRODUCTION = ones(NG,1)*NaN;
            else                                                            % othervise, calculate the transition cost
                STATE_DIFF = CURRENT_STATE - PREVIOUS_STATE;
                % STATE_DIFF = 1  means unit is commited
                % STATE_DIFF = -1 means unit is decommited
                if START_UP_COST_METHOD == 1   % start-up costs are constant and equal to cold start cost
                    GEN_START_SHUT_COST(:,K) = (STATE_DIFF > 0) .* GSC;
                elseif START_UP_COST_METHOD == 2
                    GEN_START_SHUT_COST(:,K) =                            ((STATE_DIFF > 0) & (-X_PREV(:,K) >= (GMINDOWN + GCSTIME))) .* GSC;  % cold start-up cost
                    GEN_START_SHUT_COST(:,K) = GEN_START_SHUT_COST(:,K) + ((STATE_DIFF > 0) & (-X_PREV(:,K) <  (GMINDOWN + GCSTIME))) .* GSH;  % hot start-up cost
                else
                    GEN_START_SHUT_COST(:,K) = (STATE_DIFF > 0) .* (ALPHA + BETA .* (1-exp(X_PREV(:,K) ./ TAU)));   % exponential start-up costs
                end
                GEN_START_SHUT_COST(:,K) = GEN_START_SHUT_COST(:,K) + (STATE_DIFF  < 0) .* GSDC;   % shut down cost

                % find the generation [MW] and production cost for each unit
                [GEN_PRODUCTION,PROD_COST] = production(CURRENT_STATE,PREVIOUS_STATE,GMIN,GMAX,DEMAND,HOUR,GNLC,GFC,GINC,NG,GRAMPUP,GRAMPDOWN,PRODUCTION_PREV(:,K),GEN_ORDER,COEF_A,COEF_B,COEF_C,DISPATCH_METHOD);
            end
            X_TEMP(:,K) = X; % save the updated gen. work. times when moved from previous state to the current one
            PRODUCTION_TEMP(:,K) = GEN_PRODUCTION; % also save the updated gen. outputs
            if HOUR == 1
                TOTAL_COST(K) = sum(PROD_COST) + sum(GEN_START_SHUT_COST(:,K));
            else
                TOTAL_COST(K) = sum(PROD_COST) + sum(GEN_START_SHUT_COST(:,K)) + FC_PREV(K);
            end % if HOUR
        end  % K

        % among all transitions from each feasible state at previous hour
        % to the current state (at current hour), save up to MN with minimal total cost
        [MM,II] = sort(TOTAL_COST(TOTAL_COST ~= 0));
        for K = 1:MN
            if MM(K) ~= Inf
                COUNTER = COUNTER +1;
                FC(COUNTER,1) = MM(K);
                TR(COUNTER,1:size(TR_PREV,2)) = TR_PREV(II(K),:);
                TR(COUNTER,end) = FEASIBLE_STATES_NUM(J);
                X_CURR(:,COUNTER) = X_TEMP(:,II(K));
                PRODUCTION_CURR(:,COUNTER) = PRODUCTION_TEMP(:,II(K));
            end % if MM(K)
        end % if K
    end   % J

    if COUNTER == 0;                                                                    % If the rest of the list is empty, then it means
        STR = ['NO FEASIBLE STATES FOR HOUR ',num2str(HOUR),'! PROGRAM TERMINATES!'];   % there are no feasible states,
        msgbox(STR,'NO FEASIBLE STATES','warn');                                        % and program terminates
        return
    end
end   % HOUR

%============================================
% END OF SEARCHING PROCEDURE
% ============================================
% The search is complete. Now program finds the best solution (least expensive state) at the last hour of the optimization horizon.
[M,I]=min(FC(1:COUNTER));
BEST_PATH = TR(I,:).';    % find the best transition path
% evaluate the solution and print the results
evaluate_solution(NT,BEST_PATH,LIST_STATES,GMIN,GMAX,DEMAND,GEN_ORDER,GNLC,GFC,GINC,GSC,INI_STATE,NG,...
    GRAMPUP,GRAMPDOWN,COEF_A,COEF_B,COEF_C,DISPATCH_METHOD,DETAIL_PRINT_FLAG,GSDC,GSTATINI,...
    GMINUP,GMINDOWN,START_UP_COST_METHOD,GCSTIME,GSH,ALPHA,BETA,TAU)
warning on
t=toc;
fprintf('\n Elapsed time: %15.4f sec.\n\n',t)
end

%===============================================================================================================
% FUNCTIONS DEFINITIONS
%===============================================================================================================
function [GMAXcum,GMINcum,LIST_STATES,LIST_INDEX] = priority_list(GNLC,GFC,GMAX,GMIN,GINC,NG)
%% --------------------------------------------------------------------------------------------
% Creates the list of states, where each state has one unit commited more than
% the previous state. Generators are ordered according to their full load average cost,
% the least expensive coming first.
% The list is column-based: the 1st column contains only the cheapest gen.;
% the 2nd column contains two cheapest etc.; the last one containts all NG gen.
% OUTPUT:
% LIST_STATES [NG x NG] - matrix of states; each column represents one state
% LIST_INDEX [NG x 1]   - order of generator indices (least expensive gen. is the first)
% GMINcum [NG x 1]      - total min. generator output for each state
% GMAXcum [NG x 1]      - total max. generator output for each state
% Example: if the order of generators average costs gives G3,G1,G2, then:
% LIST_INDEX = 3,1,2
% LIST_STATES = [0 1 1       - 1st column has only          G3 commited
%                0 0 1       - 2nd column has commited      G3+G1
%                1 1 1]      - 3rd column has commited      G3+G1+G2
%--------------------------------------------------------------------------------------------
GFULLAVECOST = (0*GNLC + GFC.*GMAX.*GINC/1000)./GMAX;   % Calculate full load average cost for each unit
[M,LIST_INDEX] = sort(GFULLAVECOST);                    % sort them (make a priority list)
LIST_STATES = triu(ones(NG));                           % prepare matrix of the states [NGxNG]
LIST_STATES(LIST_INDEX,:) = LIST_STATES(1:NG,:);        % create the list of states
LIST_STATES = logical(LIST_STATES);                     % not neccessarily, but it is a good manner
GMAXcum = cumsum(GMAX(LIST_INDEX));                     % create total max. generator output for each state
GMINcum = cumsum(GMIN(LIST_INDEX));                     % create total min. generator output for each state

prints_states(NG,GMINcum,GMAXcum,LIST_STATES)
end

function [GMAXlst,GMINlst,LIST_STATES,GEN_ORDER] = complete_enumeration(GNLC,GFC,GMAX,GMIN,GINC,NG)
%% --------------------------------------------------------------------------------------------
% Creates the complete list of states (for NG generators, there are totally 2^NG states)
% The list is column-based (each column represents a state).
% OUTPUT:
% LIST_STATES [NG x 2^NG] - matrix of states; each column represents one state
% GEN_ORDER [NG x 1]      - order of generators to be commited (least expensive gen. the first)
% GMINlst [2^NG x 1]      - total min. generator output for each state
% GMAXlst [2^NG x 1]      - total max. generator output for each state
% Example: if there are 3 generators, then there are totally 8 states:
% LIST_STATES = [0 0 0 0 1 1 1 1
%                0 0 1 1 0 0 1 1
%                0 1 0 1 0 1 0 1]
%--------------------------------------------------------------------------------------------
GFULLAVECOST = (0*GNLC + GFC.*GMAX.*GINC/1000)./GMAX; % Calculate full load average cost for each unit
[M,GEN_ORDER] = sort(GFULLAVECOST);                   % sort them (make a priority list of gen. commitment)
LIST_STATES=dec2bin(0:2^NG-1)';                       % all possible combinations of NG generators
LIST_STATES = logical(sscanf(LIST_STATES,'%1d',size(LIST_STATES)));

GMINlst = LIST_STATES.' * GMIN; % for each state (combination of generators) in the list,
GMAXlst = LIST_STATES.' * GMAX; % find the min. and max. power of the combination

% next 3 lines are not neccessary, but it is nice to see max. posisble output of generators in increasing order
[GMAXlst,INDEX]=sort(GMAXlst);      % sort the states according to increasing total max. power
GMINlst = GMINlst(INDEX);           % re-order min. power accordingly
LIST_STATES = LIST_STATES(:,INDEX); % and re-order the list of states as well

prints_states(NG,GMINlst,GMAXlst,LIST_STATES)
end

function [FEASIBLE_STATES_NUM,SUCCESS] = find_feasible_states(GMINlst,GMAXlst,DEMAND,HOUR,RES_UP,RES_DN)
%% --------------------------------------------------------------------------------------------------------------
% Determines all feasible states from the list of possible states
% Feasable states are the states where demand is between total min and total max output of commited generators
% If no feasible states found, program prepares termination
% OUTPUT:
% FEASIBLE_STATES_NUM   - vector of positions (columns) of feasible states in the list of states for current hour
% SUCCESS               - indicator: 1 - found at least one feasible states; 0 - no feasible states found
%----------------------------------------------------------------------------------------------------------------
FEASIBLE_STATES_NUM = find((GMINlst <= DEMAND(HOUR)-RES_DN(HOUR)) & (DEMAND(HOUR)+RES_UP(HOUR) <= GMAXlst));

if isempty(FEASIBLE_STATES_NUM)         % if there are no feasible states found
    SUCCESS = 0;                        % prepare for program termination
    STR = ['NO FEASIBLE STATES FOR HOUR ',num2str(HOUR),'! PROGRAM TERMINATES!'];
    msgbox(STR,'NO FEASIBLE STATES','warn');
    return
else
    SUCCESS = 1;
end
end

function [GENERATION,PROD_COST] = production(CURRENT_STATE,PREVIOUS_STATE,GMIN,GMAX,DEMAND,HOUR,GNLC,GFC,GINC,NG,GRAMPUP,GRAMPDOWN,PRODUCTION_PREV,GEN_ORDER,COEF_A,COEF_B,COEF_C,DISPATCH_METHOD)
%% --------------------------------------------------------------------------------------------------------------
% For the given state, calculates the MW output for each commited generator
% so the total production costs are minimal.
% INPUT: DISPATCH_METHOD
% DISPATCH_METHOD = 3 - uses quick linear generator dispatch (one segment cost curve)
% DISPATCH_METHOD = 2 - same as previous, just using linprog from optimization toolbox
% DISPATCH_METHOD = 1 - generator dispatch using quadprog from optimization toolbox;
%                          in this case cost curve is quadratic: GEN_COST = A + B*PG + C*PG^2
%                          where A,B,C are the coefficents and PG is generator output.
% OUTPUT:
% GENERATION [NG x 1] - vector of power output for each generator
% PROD_COST  [NG x 1] - vector of generation cost for each generator
% ---------------------------------------------------------------------------------------------------------------
if DISPATCH_METHOD == 3                          % DISPATCH_METHOD = 3 means quick linear dispatch
    lb = zeros(size(GMIN));
    ub = zeros(size(GMAX));
    if HOUR ==1
        lb = GMIN .* CURRENT_STATE;                 % lower bounds for generator output
        ub = GMAX .* CURRENT_STATE;                 % upper bounds for generator output
    else
        lb(CURRENT_STATE == 1) = max([GMIN(CURRENT_STATE == 1),PRODUCTION_PREV(CURRENT_STATE == 1)-GRAMPDOWN(CURRENT_STATE == 1)],[],2)...
            .* CURRENT_STATE(CURRENT_STATE == 1);                 % upper bounds for generator output

        ub(PREVIOUS_STATE == 0) = min([GMAX(PREVIOUS_STATE == 0),max([GRAMPUP(PREVIOUS_STATE == 0),GMIN(PREVIOUS_STATE == 0)],[],2)],[],2)...
            .* CURRENT_STATE(PREVIOUS_STATE == 0);
        ub(PREVIOUS_STATE == 1) = min([GMAX(PREVIOUS_STATE == 1),PRODUCTION_PREV(PREVIOUS_STATE == 1)+GRAMPUP(PREVIOUS_STATE == 1)],[],2)...
            .* CURRENT_STATE(PREVIOUS_STATE == 1);                 % upper bounds for generator output
    end

    if (sum(lb) > DEMAND(HOUR)) | (sum(ub) < DEMAND(HOUR)) | any((ub-lb) < 0)
        GENERATION = ones(NG,1)*NaN;
        PROD_COST = ones(NG,1)*Inf;
    else
        GENERATION = dispatch(CURRENT_STATE,lb,ub,DEMAND,HOUR,GEN_ORDER);
        PROD_COST =  GNLC .* CURRENT_STATE + GFC .* GINC .* GENERATION .* CURRENT_STATE / 1000; % and calculate their costs
    end
    return
else
    Aeq = double(CURRENT_STATE.');              % sum of output of commited generators
    beq = DEMAND(HOUR);                         % must match demand

    lb = zeros(size(GMIN));
    ub = zeros(size(GMAX));
    if HOUR ==1
        lb = GMIN .* CURRENT_STATE;                 % lower bounds for generator output
        ub = GMAX .* CURRENT_STATE;                 % upper bounds for generator output
    else
        lb(CURRENT_STATE == 1) = max([GMIN(CURRENT_STATE == 1),PRODUCTION_PREV(CURRENT_STATE == 1)-GRAMPDOWN(CURRENT_STATE == 1)],[],2)...
            .* CURRENT_STATE(CURRENT_STATE == 1);                 % upper bounds for generator output

        ub(PREVIOUS_STATE == 0) = min([GMAX(PREVIOUS_STATE == 0),max([GRAMPUP(PREVIOUS_STATE == 0),GMIN(PREVIOUS_STATE == 0)],[],2)],[],2)...
            .* CURRENT_STATE(PREVIOUS_STATE == 0);
        ub(PREVIOUS_STATE == 1) = min([GMAX(PREVIOUS_STATE == 1),PRODUCTION_PREV(PREVIOUS_STATE == 1)+GRAMPUP(PREVIOUS_STATE == 1)],[],2)...
            .* CURRENT_STATE(PREVIOUS_STATE == 1);                 % upper bounds for generator output
    end
    options = optimset('Display','Off');        % supress displays of linprog function
    if DISPATCH_METHOD == 2                      % economic dispatch using linprog
        f = GFC .* GINC .* CURRENT_STATE / 1000;    % vector of fuel cost for each generator
        [GENERATION,FVAL,EXITFLAG] = linprog(f,[],[],Aeq,beq,lb,ub,[],options);      % calculate the optimal production for each generator
        if EXITFLAG > 0
            PROD_COST =  GNLC .* CURRENT_STATE + GFC .* GINC .* GENERATION .* CURRENT_STATE / 1000; % and calculate their costs
        else
            GENERATION = ones(NG,1)*NaN;
            PROD_COST = ones(NG,1)*Inf;
        end
    end
    if DISPATCH_METHOD == 1                      % economic dispatch using quadrog
        %-------------------------------------------------------------------------
        % If quadprog is called with all set of generators, no matter if they are commited or not,
        % then the dimension of the quadratic programming problem is NG (NG decision variables), even
        % though some of generators can not generate (lb = ub = 0). This approach can lead quadprog to
        % infeasible solution (EXITFLAG = -2). With some experimenting, it is concluded that initial
        % guess is of crucial importance for this case. If X0 is determined by a quick dispatch function,
        % no infeasibility has reported.
        %-------------------------------------------------------------------------
        %         H = 2*diag(COEF_C .* CURRENT_STATE);
        %         f = COEF_B .* CURRENT_STATE;    % vector of fuel cost for each generator
        %         X0 = dispatch(CURRENT_STATE,GMIN,GMAX,DEMAND,HOUR,GEN_ORDER);                   % find approximate initial conditions
        %         [GENERATION,FVAL,EXITFLAG] = quadprog(H,f,[],[],Aeq,beq,lb,ub,X0,options);      % calculate the optimal production for each generator
        %         if EXITFLAG > 0
        %             PROD_COST =  COEF_A.*CURRENT_STATE + COEF_B.*GENERATION.*CURRENT_STATE + COEF_C.*GENERATION.^2.*CURRENT_STATE; % and calculate their costs
        %         else
        %             GENERATION = ones(NG,1)*NaN;
        %             PROD_COST = ones(NG,1)*Inf;
        %         end
        %-------------------------------------------------------------------------
        % However, in order to avoid call to a quick dispatch function and therefore to speed up
        % the program, only commited generators should be provided to quadprog. (where CURRENT_STATE=1)
        % In this case, the size of the problem reduces (size <= NG) and no infeasibility has encountered.
        GENERATION = zeros(NG,1);
        X0 = [];
        H = 2*diag(COEF_C(CURRENT_STATE));
        f = COEF_B(CURRENT_STATE);
        Aeq = Aeq(:,CURRENT_STATE);
        lb = lb(CURRENT_STATE);
        ub = ub(CURRENT_STATE);
        [GENERATION1,FVAL,EXITFLAG] = quadprog(H,f,[],[],Aeq,beq,lb,ub,X0,options);      % calculate the optimal production for each generator
        if EXITFLAG > 0
            GENERATION(CURRENT_STATE) = GENERATION1;
            PROD_COST =  (COEF_A.*CURRENT_STATE) + (COEF_B.*GENERATION.*CURRENT_STATE) + (COEF_C.*GENERATION.^2.*CURRENT_STATE); % and calculate their costs
        else
            GENERATION = ones(NG,1)*NaN;
            PROD_COST = ones(NG,1)*Inf;
        end
        %-------------------------------------------------------------------------
    end
end
end

function prints_states(NG,GMINcum,GMAXcum,LIST_STATES)
%% --------------------------------------------------------------------------------------------------------------
% Prints out the list of all possible states
% Note that priority list consists of NG states and the enumeration lists contains 2^NG states
% ---------------------------------------------------------------------------------------------------------------
fprintf('   State No.      MW min        MW max                     Units\n')
fprintf('%s',repmat(' ',1,23))
fprintf(['               ',repmat('    %5d ', 1, NG)],1:NG)
fprintf('\n %s \n',repmat('-',1,80'))
for I=1:size(LIST_STATES,2)
    fprintf('      %2d       %8.1f      %8.1f ',I,GMINcum(I),GMAXcum(I))
    fprintf([repmat('       %2d ', 1, size(LIST_STATES,1)) '\n'], LIST_STATES(:,I));
end
end


function evaluate_solution(NT,BEST_PATH,LIST_STATES,GMIN,GMAX,DEMAND,GEN_ORDER,GNLC,GFC,GINC,GSC,INI_STATE,NG,...
    GRAMPUP,GRAMPDOWN,COEF_A,COEF_B,COEF_C,DISPATCH_METHOD,DETAIL_PRINT_FLAG,GSDC,GSTATINI,...
    GMINUP,GMINDOWN,START_UP_COST_METHOD,GCSTIME,GSH,ALPHA,BETA,TAU)
%% --------------------------------------------------------------------------------------------------------------
% For the given set BEST_PATH of the states in each time step, calculates production costs, transition costs
% total costs etc.
% In the end, this function calls function to print the results in a tabulated form.
% This function follows pretty much the same fashion as the main procedure and will not be described in details.
% The only difference is that now the optimal path is known, so there is no need for searching.
% For each state in the path this procedure calculates the costs and finally calls the printing routine.
% ---------------------------------------------------------------------------------------------------------------
GEN_START_SHUT_COST1 = zeros(NG,NT);
GEN_PRODUCTION1      = zeros(NG,NT);
PROD_COST1           = zeros(NG,NT);
FCOST1               = zeros(NT,1);
GENERATING_COST1     = zeros(NT,1);
GEN_PRODUCTION       = zeros(NG,1);
X  = GSTATINI;
for HOUR = 1:NT
    PREV_STATES_NUM = BEST_PATH(HOUR);
    FEASIBLE_STATES_NUM = BEST_PATH(HOUR+1);
    X_PREV = X;
    if HOUR==1 & PREV_STATES_NUM == 0
        PREVIOUS_STATE = INI_STATE;
    else
        PREVIOUS_STATE = LIST_STATES(:,PREV_STATES_NUM);
    end
    CURRENT_STATE  = LIST_STATES(:,FEASIBLE_STATES_NUM);
    PRODUCTION_PREV = GEN_PRODUCTION;
    [GEN_PRODUCTION,PROD_COST] = production(CURRENT_STATE,PREVIOUS_STATE,GMIN,GMAX,DEMAND,HOUR,GNLC,GFC,GINC,NG,GRAMPUP,GRAMPDOWN,PRODUCTION_PREV,GEN_ORDER,COEF_A,COEF_B,COEF_C,DISPATCH_METHOD);

    STATE_DIFF = CURRENT_STATE - PREVIOUS_STATE;
    [X,SUCCESS] = check_up_down_time(CURRENT_STATE,PREVIOUS_STATE,X_PREV,GMINUP,GMINDOWN,NG);


    if START_UP_COST_METHOD == 1   % start-up costs are constant and equal to cold start cost
        GEN_START_SHUT_COST = (STATE_DIFF > 0) .* GSC;
    elseif START_UP_COST_METHOD == 2
        GEN_START_SHUT_COST =                       ((STATE_DIFF > 0) & (-X_PREV >= (GMINDOWN + GCSTIME))) .* GSC;  % hot start-up cost
        GEN_START_SHUT_COST = GEN_START_SHUT_COST + ((STATE_DIFF > 0) & (-X_PREV <  (GMINDOWN + GCSTIME))) .* GSH;  % cold start-up cost
    else
        GEN_START_SHUT_COST = (STATE_DIFF > 0) .* (ALPHA + BETA .* (1-exp(X_PREV ./ TAU)));
    end

    GEN_START_SHUT_COST = GEN_START_SHUT_COST + (STATE_DIFF < 0 ) .* GSDC;       % shut down cost

    if HOUR == 1
        TOTAL_COST = sum(PROD_COST) + sum(GEN_START_SHUT_COST);
    else
        TOTAL_COST = sum(PROD_COST) + sum(GEN_START_SHUT_COST) + FCOST1(HOUR-1);
    end % if HOUR
    FCOST1(HOUR) = TOTAL_COST;
    GENERATING_COST1(HOUR) = sum(PROD_COST);
    GEN_PRODUCTION1(:,HOUR) = GEN_PRODUCTION;
    PROD_COST1(:,HOUR) = PROD_COST;
    GEN_START_SHUT_COST1(:,HOUR) = GEN_START_SHUT_COST;

end   % HOUR = 1:NT
GEN_START_SHUT_COST_TOTAL = sum(GEN_START_SHUT_COST1).';

print_results(BEST_PATH,LIST_STATES,INI_STATE,NT,NG,GMIN,GMAX,DEMAND,FCOST1,GENERATING_COST1,GEN_PRODUCTION1,PROD_COST1,GEN_START_SHUT_COST1,DETAIL_PRINT_FLAG)
end

function print_results(BEST_PATH,LIST_STATES,INI_STATE,NT,NG,GMIN,GMAX,DEMAND,FCOST1,GENERATING_COST1,GENERATION,PROD_COST,START_COST,DETAIL_PRINT_FLAG)
if DETAIL_PRINT_FLAG == 0
    S = ['Hour         '
        'Demand       '
        'Tot.Gen      '
        'Min MW       '
        'Max MW       '
        'ST-UP Cost   '
        'Prod.Cost    '
        'F-Cost       '
        'State        '
        'Units ON/OFF '];
    fprintf('\n%s',repmat('=',1,150'))
    fprintf('\n       HOURLY RESULTS:')
    fprintf('\n%s \n',repmat('=',1,150'))
    fprintf([repmat('%12s ', 1, size(S,1))], S');
    fprintf('\n%s\n',repmat('-',1,150'))
else
    S = ['UNITS          '
        'ON/OFF         '
        'GENERATION     '
        'MIN MW         '
        'MAX MW         '
        'ST-UP Cost     '
        'PROD.COST      '];
end

if BEST_PATH(1) == 0
    LIST_STATES = [LIST_STATES,INI_STATE];
    BEST_PATH(1) = size(LIST_STATES,2);
end

for HOUR = 1:length(BEST_PATH)-1
    CURRENT_STATES_NUM  = BEST_PATH(HOUR+1);
    CURRENT_STATE   = LIST_STATES(:,CURRENT_STATES_NUM);

    MIN_MW = CURRENT_STATE.*GMIN;
    MAX_MW = CURRENT_STATE .*GMAX;
    if HOUR ==1 & DETAIL_PRINT_FLAG == 0
        fprintf('%3d  %12s  %12s %12.0f %12.0f %12.0f  %12.0f %12.0f %10.0f ',HOUR-1, '-','-',sum(LIST_STATES(:,BEST_PATH(HOUR)).*GMIN),sum(LIST_STATES(:,BEST_PATH(HOUR)).*GMAX),0,0,0,BEST_PATH(HOUR))
        fprintf(['       ',repmat('%2d', 1, size(LIST_STATES(:,BEST_PATH(HOUR)),1)),'\n'], LIST_STATES(:,BEST_PATH(HOUR)));
    end

    if DETAIL_PRINT_FLAG == 0;
        fprintf('%3d  %12.0f  %12.0f %12.0f %12.0f ',HOUR, DEMAND(HOUR), sum(GENERATION(:,HOUR)), sum(MIN_MW), sum(MAX_MW))
        fprintf('%12.0f  %12.0f %12.0f %10d ',sum(START_COST(:,HOUR)),sum(PROD_COST(:,HOUR)),FCOST1(HOUR),CURRENT_STATES_NUM)
        fprintf(['       ',repmat('%2d', 1, size(CURRENT_STATE,1)),'\n'], CURRENT_STATE);
    else
        TEMP = [(1:NG).',CURRENT_STATE,GENERATION(:,HOUR),MIN_MW,MAX_MW,START_COST(:,HOUR),PROD_COST(:,HOUR)];
        fprintf('\n\n\nHOUR: %2d             DEMAND:%7.1f MW           F-COST: %6.1f £',HOUR,DEMAND(HOUR),FCOST1(HOUR))
        fprintf('\n%s \n',repmat('-',1,120'))
        fprintf([repmat('%15s ', 1, size(S,1)) '\n\n'], S');fprintf('\n');
        fprintf(['%3d %15d ',repmat('%15.1f', 1, size(TEMP,2)-2) '\n'], TEMP.');
        fprintf('%s \n',repmat('-',1,120'))
        fprintf('TOTAL: %12d  %14.1f %14.1f %14.1f %14.1f %14.1f\n',sum(CURRENT_STATE),sum(GENERATION(:,HOUR)), sum(MIN_MW), sum(MAX_MW),sum(START_COST(:,HOUR)),sum(PROD_COST(:,HOUR)))
    end
end
end

function [X_CURR,SUCCESS] = check_up_down_time(CURRENT_STATE,PREVIOUS_STATE,X_PREV,GMINUP,GMINDOWN,NG)
%% --------------------------------------------------------------------------------------------------------------
% Checks wether the transition from previous state to the current state is feasible
% from the minimum up and down times point of view.
% OUTPUT:
% X_CURR [NG x 1]   - vector of working hours for the new state (NaN if transition is not possible)
% SUCCESS           - indicator: 1 - transition is possible; 0 - transition is not possible
%----------------------------------------------------------------------------------------------------------------
X_CURR = zeros(NG,1 );
SUCCESS = 1;
% for the current state of generators, first check if any generator
% has been ON less than GMINUP or been OFF less than GMINDOWN
if all((X_PREV - GMINUP).*(PREVIOUS_STATE - CURRENT_STATE) >=0 & (-X_PREV - GMINDOWN).*(CURRENT_STATE - PREVIOUS_STATE) >=0)
    for I=1:NG
        % current state is feasible regarding min up and down times; now calculate X_CURR - working times for each unit
        if (X_PREV(I) >= 1) & (CURRENT_STATE(I) == 1)
            X_CURR(I) = X_PREV(I) + 1;
        elseif (X_PREV(I) <= -GMINDOWN(I)) & (CURRENT_STATE(I) == 1)
            X_CURR(I) = 1;
        elseif (X_PREV(I) <= -1) & (CURRENT_STATE(I) == 0)
            X_CURR(I) = X_PREV(I) - 1;
        elseif (X_PREV(I) >= GMINUP(I)) & (CURRENT_STATE(I) == 0)
            X_CURR(I) = -1;
        end
    end
else                                % current state violates min up and down times
    SUCCESS = 0;                    % set the indicator to zero (failed),
    X_CURR = ones(NG,1 )*NaN;       % also set the working times to NaNs
    return                          % and stop further working time calculation
end
end

function GENERATION = dispatch(CURRENT_STATE,GMIN,GMAX,DEMAND,HOUR,GEN_ORDER)
%% --------------------------------------------------------------------------------------------------------------
% For the given state, calculates the MW output for each commited generator
% Generators are dispatched in a merit order (first the least expensive, last the most expensive)
% Note: GEN_ORDER is based on No Load Cost, Fuel Cost and Incremental costs.
% OUTPUT:
% GENERATION [NG x 1] - vector of power output for each generator
% ---------------------------------------------------------------------------------------------------------------
GENERATION = GMIN.*CURRENT_STATE;               % first set the output for each commited generator to their minimal stable generation
LOAD = DEMAND(HOUR) - sum(GENERATION);          % then reduce the load for the total minimal generation
for K = 1:length(CURRENT_STATE);                % note that CURRENT_STATE is the feasible one, ie.  demand may be supplied by commited generators
    L = GEN_ORDER(K);                                                           % GEN_ORDER is the merit list of dispatching generators
    GENERATION(L) = GENERATION(L) + min(GMAX(L)-GMIN(L),LOAD)*CURRENT_STATE(L); % increase the power of the next generator in the list
    LOAD = LOAD - min(GMAX(L)-GMIN(L),LOAD)*CURRENT_STATE(L);                   % either to their max. or to match the load
end                                                                             % whenever generation increases, load reduces, until they match
end