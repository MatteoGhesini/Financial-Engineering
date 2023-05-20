function Display_Result(AO,RC,VR,EL)
%Common function to the four method to display all results
%
%

disp('Add-On: ')
disp(['LGD     (only):  ', num2str(round(AO(1)*100,2)), '% | ', num2str(round(AO(2)*100,2)), '%'])
disp(['k       (only):  ', num2str(round(AO(3)*100,2)), '% | ', num2str(round(AO(4)*100,2)), '%'])
disp(['LGD,k   (ind):   ', num2str(round(AO(5)*100,2)), '% | ', num2str(round(AO(6)*100,2)), '%'])
disp(['LGD,k   (cor):   ', num2str(round(AO(7)*100,2)), '% | ', num2str(round(AO(8)*100,2)), '%'])
disp(' ')
disp('RC: ')
disp(['LGD     (only):  ', num2str(round(RC(1),4)), ' | ', num2str(round(RC(2),4)) ])
disp(['k       (only):  ', num2str(round(RC(3),4)), ' | ', num2str(round(RC(4),4)) ])
disp(['LGD,k   (ind):   ', num2str(round(RC(5),4)), ' | ', num2str(round(RC(6),4)) ])
disp(['LGD,k   (cor):   ', num2str(round(RC(7),4)), ' | ', num2str(round(RC(8),4)) ])
disp(' ')
disp('VaR: ')
disp(['LGD     (only):  ', num2str(round(VR(1),4)), ' | ', num2str(round(VR(2),4)) ])
disp(['k       (only):  ', num2str(round(VR(3),4)), ' | ', num2str(round(VR(4),4)) ])
disp(['LGD,k   (ind):   ', num2str(round(VR(5),4)), ' | ', num2str(round(VR(6),4)) ])
disp(['LGD,k   (cor):   ', num2str(round(VR(7),4)), ' | ', num2str(round(VR(8),4)) ])
disp(' ')
disp('EL: ')
disp(['LGD     (only):  ', num2str(round(EL(1),4)), ' | ', num2str(round(EL(2),4)) ])
disp(['k       (only):  ', num2str(round(EL(3),4)), ' | ', num2str(round(EL(4),4)) ])
disp(['LGD,k   (ind):   ', num2str(round(EL(5),4)), ' | ', num2str(round(EL(6),4)) ])
disp(['LGD,k   (cor):   ', num2str(round(EL(7),4)), ' | ', num2str(round(EL(8),4)) ])

end