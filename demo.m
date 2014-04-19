% Create some test matrices
R1_dallen = [...
    -0.4129     0.8744      0.2547  ;...
	-0.6783     -0.1086     -0.7267 ;...
	-0.6077     -0.4729     0.6380  ...
    ]
    
R2_dallen = [...
    -0.6101     0.5428      0.5772  ;...
    -0.4169     0.3996      -0.8164 ;...
    -0.6738     -0.7387     -0.0175 ...
    ]

w_dallen = 1.0/3.0

Rslerped_dallen = slerpRotationMatrices(R1_dallen, R2_dallen, w_dallen);