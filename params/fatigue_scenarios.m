function S = fatigue_scenarios(name)
% FATIGUE_SCENARIOS  Fatigue scenario parameters for Singh & Sharma, IEEE CSL 2025.
%
%   S = fatigue_scenarios(NAME) returns a struct S for scenario NAME.
%
%   NAME options:
%       'A_mild'        -- Case A: Mild gradual fatigue (maintenance therapy)
%       'B_rapid'       -- Case B: Rapid exhaustion    (intensive therapy)
%       'C_prefatigued' -- Case C: Pre-fatigued start  (residual fatigue)
%
%   Fields returned:
%       S.phi0      -- Initial FES active fraction  [1 x nf]
%       S.wf_mult   -- Fatigue rate multiplier  (applied to wf0)
%       S.wr_mult   -- Recovery rate multiplier (applied to wr0)
%       S.label     -- Short display label
%       S.description -- Clinical description

switch lower(name)

    case 'a_mild'
        S.phi0      = ones(1, 5);      % Fully rested
        S.wf_mult   = 0.7;
        S.wr_mult   = 1.3;
        S.label     = 'A – Mild gradual';
        S.description = 'Slow decline; maintenance therapy with rested patient.';

    case 'b_rapid'
        S.phi0      = ones(1, 5);      % Fully rested at start
        S.wf_mult   = 2.5;
        S.wr_mult   = 0.4;
        S.label     = 'B – Rapid exhaustion';
        S.description = 'Fast depletion; intensive therapy driving toward near-complete fatigue.';

    case 'c_prefatigued'
        S.phi0      = 0.45 * ones(1, 5);  % Patient arrives pre-fatigued
        S.wf_mult   = 1.5;
        S.wr_mult   = 0.8;
        S.label     = 'C – Pre-fatigued start';
        S.description = 'Residual fatigue at session start; most challenging clinical scenario.';

    otherwise
        error('fatigue_scenarios: unknown scenario ''%s''. Use A_mild, B_rapid, or C_prefatigued.', name);
end

end
