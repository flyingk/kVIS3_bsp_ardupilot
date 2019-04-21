plot_definition = {
... Row, Col, AxesLayout, PlotStyle, Color,     Group,        Channel, ScaleFactor, UnitOverride;
      1,   1,        'S',        '',    '',    'Base',           'Az',            1,          '';
      2,   1,        'S',        '',    '',    'Base',           'Gy',       180/pi,     'deg/s';
      3,   1,        'S',        '',    '',    'Base',          'AoA',            1,          '';
      4,   1,        'S',        '',    '',    'Base',        'Pitch',            1,          '';
      5,   1,        'L',        '',    '',    'Base',     'Elevator',            1,          '';
      5,   1,        'R',        '',    '', 'Flight Control 1', 'Ch3 CMD (el)',            1,          '';
    };
CUSTOM_PLOTS.(BSP_NAME).Base = plot_definition;
