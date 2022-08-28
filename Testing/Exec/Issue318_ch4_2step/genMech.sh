cd ${PELE_PHYSICS_HOME}/Support/ceptr
mechFolder=${PELE_PHYSICS_HOME}/Support/Mechanism/Models/chem-CH4-2step

poetry run ck2yaml --input=${mechFolder}/mechanism.inp --thermo=${mechFolder}/therm.dat --transport=${mechFolder}/tran.dat --permissive
poetry run convert -f ${mechFolder}/mechanism.yaml
cd /Users/mhassana/Desktop/IssuePP/PelePhysics/Testing/Exec/Issue318_ch4_2step
