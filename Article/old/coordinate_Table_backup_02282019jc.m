SessionNAME = ['11_w10d4'; '12_w11d5'; '14_w14d2'; '14_w14d5'; '14_w14d7'; '14_w14d8'; '15_w08d5'; '15_w10d3'; '15_w10d7'; '15_w10d8'; '17_w10d3'; '17_w10d4'; '17_w10d5'; '17_w10d7'; '17_w10d8'  ];
Ch_Origin   = ['S1Ch8M_' ; 'S2Ch8_P' ; 'S1Ch8_P' ; 'S1Ch8_P' ; 'S1Ch8_P' ; 'S1Ch8_P' ; 'S1Ch8_P' ; 'S4Ch8MA' ; 'S1Ch8LA' ;  'S1Ch8LP'; 'S1Ch8M_' ; 'S2Ch8_L' ; 'S1Ch4_Q' ; 'S1Ch8M_' ;  'S1Ch8M_'  ];
Ch_Orient   = ['ML'      ;  'AP'     ;  'AP'     ;  'AP'     ;  'AP'     ;  'AP'     ;  'AP'     ; 'DG'      ;  'DG'     ;  'DG'     ;  'ML'     ;  'ML'     ;  'Q4'     ;    'ML'   ;   'ML'      ];
ElectrodeID = ['CC4F_n4S';'CCAE_o3S' ;'G912_o4S' ;'G912_o4S' ;'CC4F_n4S' ;'G912_o4S' ;'G912_o4S' ; 'CAED_n4S';'G912_o4S' ;'G912_o4S' ;'CAED_n4S' ; 'CCAE_o3S'; 'EOBD_o1S'; 'G912_o4S'; 'G912_o4S'  ];
AP_measured = [1.35      ;  1.58     ;  1.10     ;  1.55     ;  1.50     ;  1.55     ;  1.90     ;  1.20     ;   1.2     ;   1.5     ;  1.00     ;  1.00     ;   1.3     ;     1.45  ;     1.2     ];
ML_measured = [0.35      ;  0.90     ;  0.75     ;  0.90     ;  0.80     ;  0.85     ;  0.70     ;  0.95     ;   0.8     ;   0.7     ;  0.5      ;  1.10     ;   0.65    ;     0.65  ;     0.65    ];
DV_measured = [4.35      ;  4.30     ;  4.65     ;  4.30     ;  4.35     ;  4.90     ;  4.70     ;  4.25     ;   4.4     ;   4.6     ;  4.30     ;  4.20     ;   4.30    ;     4.30  ;     4.5     ];
APAdj       = [+0.40     ;  +0.00    ;  +0.15    ;  +0.15    ;  +0.15    ;  +0.15    ;  -0.20    ;  -0.20    ;  -0.20    ;  -0.20    ;  +0.20    ;  +0.00    ;   -0.15   ;     +0.0  ;     +0.0    ];    % due to Bending of the Headpost or electrode.
MLAdj       = [+0.10     ;  +0.00    ;  +0.00    ;  +0.00    ;  +0.00    ;  +0.00    ;  +0.20    ;  +0.15    ;  +0.15    ;  +0.25    ;  +0.15    ;  +0.15    ;   +0.15   ;     +0.15 ;    +0.15    ];    % due to Bending of the Headpost or electrode.
DVAdj       = [+0.20     ;  +0.20    ;  +0.10    ;  +0.25    ;  +0.25    ;  +0.25    ;  +0.00    ;  +0.05    ;  +0.00    ;  +0.00    ;  +0.10    ;  +0.20    ;   +0.20   ;     +0.20 ;    +0.020   ]; % due to increasing cortical lesion after several implants and some departure from the bones.

DVAdj=DVAdj-0.20; % due to difference between touch cortex and touch water. 

