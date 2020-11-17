from PyQt5 import QtCore
from PyQt5.QtCore import pyqtSlot
from PyQt5.QtWidgets import QMessageBox
import os
import pystache


# from OpenMM_Runner import OpenMMScriptRunner


class Advanced(QtCore.QThread):
    def __init__(self, parent=None):
        super(Advanced, self).__init__(parent)

    @pyqtSlot()
    def send_arg_to_Engine(self):
        rigid_water = True
        minimize = True
        no_minimize_value = True
        DCD_Reporter = True
        State_Data_Reporter = True
        Nonbounded_cutoff_active = True
        Additional_Integrator = False
        properties_active = False
        CPU_properties_active = False
        Device_ID_active = False
        precision = 'single'
        water_active = False
        # global selected_platform

        print("Simulation parameters preparing for the start ...")

        ## FORCEFIELD CONFIGURATIONS
        self.protein_ff = self.protein_forcefield_comboBox.currentText() + '.xml'

        if 'obc' or 'gbvi' in self.protein_forcefield_comboBox.currentText():
            water_active = True

        if self.protein_forcefield_comboBox.currentText() == 'charmm36':
            print("Protein FF: %s" % self.protein_ff)
            self.water_ff = '%s/%s.xml' % (self.protein_forcefield_comboBox.currentText(), self.water_forcefield_comboBox.currentText())
            print("Water FF: %s" % self.water_ff)
        else:
            self.water_ff = '%s.xml' % self.water_forcefield_comboBox.currentText()
            print("Protein FF: %s" % self.protein_ff)
            print("Water FF: %s" % self.water_ff)

        ## SOLUTION WATER MODEL
        if len((self.water_ff.split('/')[-1]).split('.')[0]) <= 5:
            self.water_model = (self.water_ff.split('/')[-1]).split('.')[0]
            print("Water Model for Solution 1: %s" % self.water_model)

        elif (self.water_ff.split('/')[-1]).split('.')[0] == 'tip4pew':
            self.water_model = (self.water_ff.split('/')[-1]).split('.')[0]
            print("Water Model for Solution 2: %s" % self.water_model)

        else:
            self.water_model = (self.water_ff.split('/')[-1]).split('.')[0][0:5]
            print("Water Model for Solution 3: %s" % self.water_model)

        try:
            self.out_dir = self.Output_Folder_textEdit.toPlainText()
            if self.out_dir == '':
                print("hiiii")
                answer = QMessageBox.question(self, 'Warning',
                                              "There is no specified Output Directory. "
                                              "If you want to specify, please click the 'Yes' button, "
                                              "otherwise the system will create an output folder named 'output'",
                                              QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
                if answer == QMessageBox.Yes:
                    return False

                if answer == QMessageBox.No:
                    from pathlib import Path
                    path_out = os.getcwd() + "/output"
                    Path(path_out).mkdir(parents=True, exist_ok=True)

                    self.out_dir = path_out
                print("HAAAAA")
        except:
            QMessageBox.critical(self, "Error", "An error occured while getting output directory")

        if not self.rigid_water_checkBox.isChecked():
            rigid_water = False

        if self.Max_minimize_iter_textEdit.toPlainText() == "":
            no_minimize_value = False

        if not self.minimize_checkBox.isChecked():
            minimize = False

        if not self.DCD_Reporter_checkBox.isChecked():
            DCD_Reporter = False

        if not self.State_Data_Reporter_checkBox.isChecked():
            State_Data_Reporter = False

        if self.nonBounded_Method_comboBox.currentText == 'NoCutoff':
            Nonbounded_cutoff_active = False

        if self.Additional_Integrators_checkBox.isChecked():
            Additional_Integrator = True

        if self.platform_comboBox.currentText() in ["CUDA", "OpenCL"]:
            properties_active = True

        if self.platform_comboBox.currentText() == "CPU":
            CPU_properties_active = True

        if self.Device_ID_checkBox.isChecked():
            Device_ID_active = True

        platform, properties, precision = Advanced_Helper_Functions.selected_platform(self, self.platform_comboBox.currentText(),
                                                                           Device_ID_active, precision)
        print("HAAAAA")
        script_structure = dict(pdb=self.upload_pdb_textEdit.toPlainText(),
                                output_folder=self.Output_Folder_textEdit.toPlainText(),

                                long_simulation_time=self.run_duration_spinBox_2.value(),
                                long_simulation_time_unit=self.long_simulation_time_unit.currentText(),

                                Number_of_CPU=self.Number_CPU_spinBox_2.value(),
                                all_cpu=self.All_CPU_checkBox.isChecked(),
                                platform=platform, properties=properties, precision=precision,
                                properties_active=properties_active, CPU_properties_active=CPU_properties_active,
                                Device_ID_active=Device_ID_active,
                                Device_Number=self.Device_Number_comboBox.currentText(),

                                protein_ff=self.protein_ff,  # FORCEFIELD
                                water_ff=self.water_ff,  # FORCEFIELD
                                model_water=self.water_model,
                                water_active=water_active,  # FORCEFIELD if protein include obc or gbvi
                                water_padding=self.water_padding_lineEdit.text(),

                                integrator_kind=self.integrator_kind_comboBox.currentText(),  # INTEGRATOR
                                integrator_time_step=self.integrator_time_step.toPlainText(),  # INTEGRATOR
                                friction=self.friction_textEdit.toPlainText(),  # ADDITIONAL INTEGRATOR
                                Temperature=self.temperature_textEdit.toPlainText(),  # ADDITIONAL INTEGRATOR
                                Additional_Integrator=Additional_Integrator,

                                NonBoundedMethod=self.nonBounded_Method_comboBox.currentText(),
                                Constraints=self.system_constraints_comboBox.currentText(),
                                Rigid_Water=rigid_water,
                                NonBounded_cutoff=self.nonbounded_CutOff_textEdit.toPlainText(),
                                Nonbounded_cutoff_active=Nonbounded_cutoff_active,

                                Number_of_steps=self.Number_of_steps_textEdit.toPlainText(),
                                Minimize=minimize,
                                no_minimize_value=no_minimize_value,
                                Max_minimization_iteration=self.Max_minimize_iter_textEdit.toPlainText(),
                                DCDReporter=DCD_Reporter,
                                State_Data_Reporter=State_Data_Reporter,
                                DCD_write_freq=self.DCD_write_freq_textEdit.toPlainText(),
                                DCD_output_name=self.DCD_Output_Name_textEdit.toPlainText(),
                                StateData_freq=self.StateData_frequency_textEdit.toPlainText(),
                                output_directory=self.out_dir
                                )
        Advanced_Helper_Functions.update_display(self, script_structure)
        # self.update_display(script_structure)


class Advanced_Helper_Functions(QtCore.QThread):
    def __init__(self, parent=None):
        super(Advanced_Helper_Functions, self).__init__(parent)

    @pyqtSlot()
    def selected_platform(self, platform_name, Device_ID_active, precision):
        print("hoooo burada")

        if platform_name == 'OpenCL' and Device_ID_active == True:
            properties = {'OpenCLPrecision': '%s' % precision,
                          'OpenCLDeviceIndex': '%s' % self.Device_Number_comboBox.currentText()}
            return platform_name, properties, precision

        if platform_name == 'OpenCL' and Device_ID_active == False:
            properties = {'OpenCLPrecision': '%s' % precision}
            return platform_name, properties, precision

        if platform_name == 'CUDA' and Device_ID_active == True:
            properties = {'CUDAPrecision': '%s' % precision, 'CUDADeviceIndex': '%s' % self.Device_Index_Number}
            return platform_name, properties, precision

        if platform_name == 'CPU':
            print("The CPU platform always uses 'mixed' precision.")
            print("Simulation process will use %s Thread(s)" % self.Number_CPU_spinBox_2.value())
            properties = {'CpuThreads': '%s' % self.Number_CPU_spinBox_2.value()}
            precision = 'mixed'
            return platform_name, properties, precision

        if platform_name == 'Reference':
            print("The Reference platform always uses 'double' precision.")
            properties = None
            precision = 'double'
            return platform_name, properties, precision
        print("System will use '%s' Platform with '%s' Precision" % (platform_name, precision))

    @pyqtSlot()
    def update_display(self, script_structure):

        # print(script_structure)
        renderer = pystache.Renderer()
        # print("ok")
        template = pystache.parse(u'''
        ###########################################################################
        ############### This script was generated by MDPerTool v0.1 ###############
        ###########################################################################
        
        from simtk.openmm import app
        from simtk.openmm.app import PME, NoCutoff, Ewald, CutoffPeriodic, CutoffNonPeriodic, HBonds, HAngles, AllBonds
        import simtk.openmm as mm
        from simtk.unit import femtosecond, picosecond, nanometer, kelvin, angstrom, atmospheres
        from sys import stdout
        from apply_pdbfixer import fix_pdb
        from simtk.openmm import *
        from mdtraj.reporters import XTCReporter

        
        print('pdb file fixing and preparing for simulation ...')
        fixed_pdb_name = fix_pdb('{{pdb}}')
        
        print('Loading pdb to simulation engine ...')
        pdb = app.PDBFile(fixed_pdb_name)

        box = pdb.topology.getUnitCellDimensions()

        print('Modeller of pdb file is preparing ...')
        modeller = mm.app.Modeller(pdb.topology, pdb.positions)
        modeller.topology.setUnitCellDimensions(box)

        print('Forcefield parameters loading to the simulation system ...')
        forcefield = app.ForceField('{{protein_ff}}'{{#water_active}}, '{{water_ff}}'{{/water_active}})
        
        print('Adding missing hydrogens to the model ...')
        modeller.addHydrogens(forcefield)

        print('Adding solvent (both water and ions) to the model to fill a rectangular box ...')
        modeller.addSolvent(forcefield, model='{{model_water}}', padding={{water_padding}})

        print('Constructing an OpenMM System')
        self.system = forcefield.createSystem(modeller.topology, nonbondedMethod={{NonBoundedMethod}},
                                              {{#Nonbounded_cutoff_active}}nonbondedCutoff={{NonBounded_cutoff}},{{/Nonbounded_cutoff_active}}
                                              constraints={{Constraints}}, rigidWater={{Rigid_Water}}, 
                                              ewaldErrorTolerance=0.005)
        
        
        self.system.addForce(mm.MonteCarloBarostat(1 * atmospheres, {{Temperature}}, 25))

        nonbonded = [f for f in self.system.getForces() if isinstance(f, NonbondedForce)][0]
        nonbonded.setUseSwitchingFunction(use=True)
        nonbonded.setSwitchingDistance(10.0 * angstrom)
        nonbonded.setUseDispersionCorrection(True)
        
        print('Creating a LangevinIntegrator.')
        integrator = mm.{{integrator_kind}}Integrator({{#Additional_Integrator}}{{Temperature}}, {{friction}}, {{/Additional_Integrator}}{{integrator_time_step}})

        platform = mm.Platform.getPlatformByName('{{platform}}')
        
        {{#properties_active}}properties = {'{{platform}}Precision': '{{precision}}'{{#Device_ID_active}},'{{platform}}DeviceIndex': '{{Device_Number}}'{{/Device_ID_active}}}{{/properties_active}}
        {{#CPU_properties_active}}properties = {'CpuThreads': '{{Number_of_CPU}}'}{{/CPU_properties_active}}
        
        simulation = app.Simulation(modeller.topology, self.system, integrator, platform{{#properties_active}}, properties{{/properties_active}})
        
        simulation.context.setPositions(modeller.positions)

        simulation.context.computeVirtualSites()
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        

        ''')

        contents = renderer.render(template, script_structure)
        print(contents)
        # OpenMMScriptRunner(script_structure)

        """


        {{#properties_active}}properties = {'{{platform}}Precision': 'mixed'{{#Device_ID_active}},'{{platform}}DeviceIndex': '{{Device_Number}}'{{/Device_ID_active}}}{{/properties_active}}


        simulation = Simulation(modeller.topology, system, integrator, platform{{#properties_active}}, properties{{/properties_active}})
        simulation.context.setPositions(modeller.positions)

        print('Minimizing...')
        {{#Minimize}}simulation.minimizeEnergy({{#no_minimize_value}}maxIterations={{/no_minimize_value}}{{Max_minimization_iteration}}){{/Minimize}}

        simulation.context.setVelocitiesToTemperature({{Temperature}})

        print('Starting')
        {{#DCDReporter}}simulation.reporters.append(DCDReporter('{{DCD_output_name}}', {{DCD_write_freq}})){{/DCDReporter}}

        {{#State_Data_Reporter}}simulation.reporters.append(StateDataReporter(stdout, {{StateData_freq}}, step=True, 
        time=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True, temperature=True, progress=True, 
        remainingTime=True, speed=True, totalSteps={{Number_of_steps}})){{/State_Data_Reporter}}

        print('Running Production...')    
        simulation.step({{Number_of_steps}})
        print('Done')
        """