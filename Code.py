############################################IMPORTS###############################################3
from PyQt5.Qt import *
from PyQt5.QtWidgets import QFileDialog, QDialog
from PyQt5 import QtCore, QtWidgets, uic
from PyQt5.QtGui import *
from numpy.lib.function_base import angle
from mainWindow import Ui_MainWindow
import pyqtgraph as pg
import vtk
from vtk.util import numpy_support
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk 
import numpy
# common packages 
import numpy as np 
import os
import copy
import matplotlib.pyplot as plt
from functools import reduce
# reading in dicom files
import pydicom
# scipy linear algebra functions 
from scipy.linalg import norm
from ipywidgets.widgets import * 
import ipywidgets as widgets
from plotly.graph_objs import *
import sys
import math as mathh
math = vtk.vtkMath()
import os
from scipy import ndimage 



################################## CTMPR CODE ##############################################

class ApplicationWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super(ApplicationWindow, self).__init__()
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        self.setGeometry(600, 300, 400, 200)
        self.setWindowTitle('Dicom Viewer')
        pg.setConfigOption('background', 'w')
        ##################Connecting UI with functions #############################################

        uis = [self.ui.original, self.ui.Coronal, self.ui.Sagittal, self.ui.Oblique]
        for ui in uis:
            ui.viewport().installEventFilter(self)
            ui.setAspectLocked()
        ############hiding the axis##################

            ui.getPlotItem().hideAxis('left')
            ui.getPlotItem().hideAxis('bottom')
            ui.setMouseEnabled(x=False, y=False)


        self.whichview=1


        self.ui.Upload_PB.clicked.connect(self.UploadImage)
        self.ui.pushButton.clicked.connect(self.show_result)
        self.ui.pushButton_2.clicked.connect(self.show_result2)
        self.ui.pushButton_3.clicked.connect(self.show_result3)
        ###########################defining the oblique angles in 3 views #################################
        self.ui.spinBox.setValue(135)
        self.ui.spinBox_2.setValue(146.5)
        self.ui.spinBox_3.setValue(146.5)
        ##########################when interacting with mouse for zoom this increments###############
        self.indexOriginal = 1
        self.indexSagittal = 1
        self.indexCoronal = 1
        self.indexOblique = 1
        ######################Arrays having the volume from different views################################


        self.ArrayDicom = []
        self.ConstPixelSpacing = []
        self.sagittalView = [] 
        self.coronalView = []
        self.obliqueView = []
        self.AxialAngle=135
        #############################Lines drawn on the views#######################################

        self.horizontalLine = pg.InfiniteLine(angle=0, movable=True, pen='r')
        self.verticalLine = pg.InfiniteLine(angle=90, movable=True, pen='b')
        self.obliqueLine = pg.InfiniteLine(angle=135, movable=True, pen='g')

        self.horizontalLineSag = pg.InfiniteLine(angle=0, movable=True, pen='r')
        self.verticalLineSag = pg.InfiniteLine(angle=90, movable=True, pen='b')
        self.obliqueLineSag = pg.InfiniteLine(angle=146.5, movable=True, pen='g')

        self.horizontalLineCor = pg.InfiniteLine(angle=0, movable=True, pen='r')
        self.verticalLineCor = pg.InfiniteLine(angle=90, movable=True, pen='b')
        self.obliqueLineCor = pg.InfiniteLine(angle=146.5, movable=True, pen='g')
        #################################when lines moves we can change in views##########################

        self.lines = [self.verticalLine, self.horizontalLine, self.verticalLineSag, self.horizontalLineSag, self.verticalLineCor, self.horizontalLineCor]
        for line in self.lines:
            line.sigPositionChanged.connect(self.Ichanged)

      

    ##############################Oblique from Axial view########################

    def show_result(self):
        """
        waits for push button, view the angle we already setted in the ui 
        then does an oblique reconstruction with the new angle
        """
        # Array2=self.ArrayDicom.copy()
        # Array2 = ndimage.rotate(Array2, angle=45, reshape=True)

        self.ui.Oblique.clear()
        self.whichview=1
        if self.ui.spinBox.value() > self.AxialAngle:
            Angletorot=self.AxialAngle-self.ui.spinBox.value()
            self.AxialAngle=self.ui.spinBox.value()
        elif self.ui.spinBox.value() < self.AxialAngle:
            Angletorot=self.AxialAngle-self.ui.spinBox.value()
            self.AxialAngle=self.ui.spinBox.value()
        else:
            Angletorot=45
            self.AxialAngle=self.ui.spinBox.value()

        self.obliqueLine.setAngle(self.AxialAngle)
        self.obliqueView=self.ObliqueReconstruct(self.ArrayDicom,Angletorot,90)
        
        for i in range(2):
            self.obliqueView=np.rot90(self.obliqueView,1,axes=(0, 1))
        self.imageOblique = pg.ImageItem(self.obliqueView[:,:,1])
        self.ui.Oblique.addItem(self.imageOblique)

        self.indexOblique=1
        self.obliqueLine.setValue((self.indexOblique*(512/self.obliqueView.shape[2]),self.indexOblique*(512/self.obliqueView.shape[2]))) # bring iso line above contrast controls
    ##############################Oblique from sagittal view########################
    def show_result2(self):
        """
        waits for push button, view the angle we already setted in the ui 
        then does an oblique reconstruction with the new angle
        """
        Angletorot=self.ui.spinBox_2.value()
        self.obliqueLineSag.setAngle(Angletorot)

        self.ui.Oblique.clear()
        self.whichview=2
        self.obliqueView=self.ObliqueReconstruct(self.ArrayDicom,180,Angletorot)
        for i in range(2):
            self.obliqueView=np.rot90(self.obliqueView,1,axes=(0, 1))
  
        self.obliqueView=self.obliqueView[:,:,::-1]
        self.obliqueView=np.flip(self.obliqueView, axis=0)
        self.imageOblique = pg.ImageItem(self.obliqueView[:,:,1])
        self.ui.Oblique.addItem(self.imageOblique)

        self.indexOblique=1
        self.obliqueLineSag.setValue((self.obliqueView.shape[2]-(self.indexOblique*(self.ArrayDicom.shape[0]/self.obliqueView.shape[2])),self.obliqueView.shape[2]-(self.indexOblique*(self.ArrayDicom.shape[2]/self.obliqueView.shape[2])))) # bring iso line above contrast controls
     ##############################Oblique from coronal view########################
    def show_result3(self):
        """
        waits for push button, view the angle we already setted in the ui 
        then does an oblique reconstruction with the new angle
        """
        Angletorot=self.ui.spinBox_3.value()
        self.obliqueLineCor.setAngle(Angletorot)
        self.ui.Oblique.clear()
        self.whichview=3
        self.obliqueView=self.ObliqueReconstruct(self.ArrayDicom,90,Angletorot)
        
        for i in range(1):
            self.obliqueView=np.rot90(self.obliqueView,1,axes=(0, 1))
        self.obliqueView=self.obliqueView[:,:,::-1]
        self.obliqueView=np.flip(self.obliqueView, axis=0)
        self.imageOblique = pg.ImageItem(self.obliqueView[:,:,1])
        self.ui.Oblique.addItem(self.imageOblique)

        self.indexOblique=1
        self.obliqueLineCor.setValue((self.obliqueView.shape[2]-(self.indexOblique*(self.ArrayDicom.shape[0]/self.obliqueView.shape[2])),self.obliqueView.shape[2]-(self.indexOblique*(self.ArrayDicom.shape[2]/self.obliqueView.shape[2])))) # bring iso line above contrast controls

##################################Function to transform from vtk to numpy array######################3
    def vtkToNumpy(self, data): 
        temp = vtk_to_numpy(data.GetPointData().GetScalars())
        dims = data.GetDimensions()
        numpy_data = temp.reshape(dims[2], dims[1], dims[0]) 
        numpy_data = numpy_data.transpose(2,1,0)    
        return numpy_data 
###################################Function to load dicom files and calculate its thickness##################
    def load_and_calculate_Thickness(self, path):
        slices = [pydicom.dcmread(path + '/' + s) for s in os.listdir(path)]
        slices = [s for s in slices if 'SliceLocation' in s]
        slices.sort(key = lambda x: int(x.InstanceNumber))
        try:
            slice_thickness = np.abs(slices[0].ImagePositionPatient[2]-slices[1].ImagePositionPatient[2])
        except:
            slice_thickness = np.abs(slices[0].SliceLocation-slices[1].SliceLocation)
        for s in slices:
            s.SliceThickness = slice_thickness
        return slices
###################takes the pixle data from metadata and construct 3d numpy array##################
    def SetImagedata(self, ArrayDicom,slices,idx):
            for s in slices:
                ArrayDicom[:,:,idx] = s.pixel_array
                idx-=1
###################################Function to know the positions of lines when dragged#######################
    def handle_sig_dragged(self,obj):
        if(obj is self.verticalLine):
            print("x = ", obj.value())
        elif(obj is self.horizontalLine):
            print("y= ", obj.value())


    def getPosAxial(self , event):
        """
        Get the position of the user click in the axial view
        """
        x = event.pos().x()
        y = event.pos().y() 

        z=self.indexOriginal

        self.verticalLine.setValue(x)
        self.horizontalLine.setValue(y)

        self.verticalLineSag.setValue(self.ArrayDicom.shape[1]-y)
        self.horizontalLineSag.setValue(self.ArrayDicom.shape[2]-z)

        self.verticalLineCor.setValue(x)
        self.horizontalLineCor.setValue(self.ArrayDicom.shape[2]-z)

    def getPosSagittal(self , event):
        """
        Get the position of the user click in the sagittal view
        """
        x = event.pos().x()      #### x Sag hena heya el y axial
        y = event.pos().y()      #### y Sag hena heya el z axial 
                                 #### z Sag hena heya el x axial shape2 of sagittal

        self.verticalLineSag.setValue(x)
        self.horizontalLineSag.setValue(y)

        z=self.indexSagittal 

        self.verticalLine.setValue(z)
        self.horizontalLine.setValue(self.sagittalView.shape[0]-x)

        self.verticalLineCor.setValue(z) 
        self.horizontalLineCor.setValue(y)

    def getPosCoronal(self , event):
        """
        Get the position of the user click in the coronal view
        """
        x = event.pos().x()     ### X coronal heya el x axial 3aks 
        y = event.pos().y()     ### Y coronal heya el z aixal 3aks
                                ### Z coronal heya el y axial 3aks
        self.verticalLineCor.setValue(x)
        self.horizontalLineCor.setValue(y)

        z=self.indexCoronal
        self.verticalLine.setValue(x)
        self.horizontalLine.setValue(z)

        self.verticalLineSag.setValue(self.coronalView.shape[2]-z)
        self.horizontalLineSag.setValue(y)

    ##############################To handle the dragging of the lines#############################

    def Ichanged(self,obj):
        """
        Cases:
        Redline Axial: Changes the slice number in the coronal view
        Blueline Axial:Changes the slice number in the sagittal view
        Redline Sagittal:Changes the slice number in axial view
        Blueline Sagittal:Changes the slice number in the coronal view
        Redline Coronal: Changes the slice number in the axial view
        Blueline Coronal: Changes the slice number in the sagittl view
        """
        if(obj is self.verticalLine):
            self.indexSagittal=int(obj.value())
            img=self.sagittalView[:,:,self.indexSagittal]  
            self.imageSagittal = pg.ImageItem(img)      
            self.ui.Sagittal.addItem(self.imageSagittal)
            self.verticalLineCor.setValue(self.indexSagittal)
        elif(obj is self.horizontalLine):
            self.indexCoronal=int(obj.value())
            img=self.coronalView[:,:,self.indexCoronal]  
            self.imageCoronal = pg.ImageItem(img)      
            self.ui.Coronal.addItem(self.imageCoronal)  
            self.verticalLineSag.setValue(self.sagittalView.shape[0]-self.indexCoronal)
        elif(obj is self.verticalLineSag):
            self.indexCoronal=self.coronalView.shape[2]-int(obj.value())
            img=self.coronalView[:,:,self.indexCoronal-1]  
            self.imageCoronal = pg.ImageItem(img)      
            self.ui.Coronal.addItem(self.imageCoronal)  
            self.horizontalLine.setValue(self.sagittalView.shape[0]-int(obj.value()))
        elif(obj is self.horizontalLineSag):
            self.indexOriginal=(self.ArrayDicom.shape[2])-int(obj.value())
            img=self.ArrayDicom[:,:,self.indexOriginal-1]  
            self.image = pg.ImageItem(img)      
            self.ui.original.addItem(self.image) 
            self.horizontalLineCor.setValue(int(obj.value()))

        elif(obj is self.verticalLineCor):
            self.indexSagittal=int(obj.value())
            img=self.sagittalView[:,:,self.indexSagittal]  
            self.imageSagittal = pg.ImageItem(img)      
            self.ui.Sagittal.addItem(self.imageSagittal)
            self.verticalLine.setValue(int(obj.value()))

        elif(obj is self.horizontalLineCor):
            self.indexOriginal=self.ArrayDicom.shape[2]-int(obj.value())
            img=self.ArrayDicom[:,:,self.indexOriginal-1]  
            self.image = pg.ImageItem(img)      
            self.ui.original.addItem(self.image) 
            self.horizontalLineSag.setValue(int(obj.value()))


    ##################Function to upload dicom Files####################
    def UploadImage(self):
        """
        The function gets the directory chosen by the user 
        and check that the directory only contains .dcm files else 
        it sends an error message and break from uploading
        then continues by constructing the colume as a numpy array to be displayed with 
        different views 
        """
        path= QFileDialog.getExistingDirectory()
        from os import listdir
        from os.path import isfile, join
        onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
        for i in onlyfiles:
            if os.path.splitext(i)[-1].lower() !='.dcm':
                error_dialog = QtWidgets.QErrorMessage()
                error_dialog.showMessage('Choose Dicom Files Only!')
                break
        slices=self.load_and_calculate_Thickness(path)



        ConstPixelDims = (int(slices[0].Rows), int(slices[0].Columns), len(slices))

        self.ConstPixelSpacing = (float(slices[0].PixelSpacing[0]), float(slices[0].PixelSpacing[1]), float(slices[0].SliceThickness))

        self.ArrayDicom = np.zeros(ConstPixelDims, dtype=slices[0].pixel_array.dtype)

        idx=len(slices)-1
        self.SetImagedata(self.ArrayDicom,slices,idx)

        for i in range(3):
            self.ArrayDicom=np.rot90(self.ArrayDicom,1,axes=(0, 1))
 
        img=self.ArrayDicom[:,:,int(self.ArrayDicom.shape[2]/2)]
                
        self.indexOriginal=int(self.ArrayDicom.shape[2]/2)

        self.image = pg.ImageItem(img)      
        self.ui.original.addItem(self.image)

        self.generateViews()
        ##################################sagittal view##############################

        self.imageSagittal = pg.ImageItem(self.sagittalView[:,:,1])
        self.ui.Sagittal.addItem(self.imageSagittal)
        
        self.ui.Sagittal.addItem(self.horizontalLineSag)
        self.ui.Sagittal.addItem(self.verticalLineSag)
        self.ui.Sagittal.addItem(self.obliqueLineSag)
        # print(mathh.degrees(mathh.atan((self.sagittalView.shape[1])/(self.sagittalView.shape[0]))))
        # self.obliqueLineSag.setAngle(90+(90-mathh.degrees(mathh.atan((self.sagittalView.shape[0])/(self.sagittalView.shape[1])))))
        ##################################coronal view##############################

        self.imageCoronal = pg.ImageItem(self.coronalView[:,:,1])
        self.ui.Coronal.addItem(self.imageCoronal)
        self.ui.Coronal.addItem(self.horizontalLineCor)
        self.ui.Coronal.addItem(self.verticalLineCor)
        self.ui.Coronal.addItem(self.obliqueLineCor)

        for i in range(2):
            self.obliqueView=np.rot90(self.obliqueView,1,axes=(0, 1))
        ##################################axial view##############################

        self.imageOblique = pg.ImageItem(self.obliqueView[:,:,0])
        self.ui.Oblique.addItem(self.imageOblique)

        self.ui.original.addItem(self.horizontalLine)
        self.ui.original.addItem(self.verticalLine)
        self.ui.original.addItem(self.obliqueLine)
        ###########################bring iso line above contrast controls##############################

        self.horizontalLine.setZValue(1000) # bring iso line above contrast controls
        self.verticalLine.setZValue(1000) # bring iso line above contrast controls
        
        self.obliqueLine.setZValue(1000) # bring iso line above contrast controls
        self.obliqueLineSag.setZValue(1000) # bring iso line above contrast controls
        self.obliqueLineCor.setZValue(1000) # bring iso line above contrast controls

        self.horizontalLineSag.setZValue(1000) # bring iso line above contrast controls
        self.verticalLineSag.setZValue(1000) # bring iso line above contrast controls

        self.horizontalLineCor.setZValue(1000) # bring iso line above contrast controls
        self.verticalLineCor.setZValue(1000) # bring iso line above contrast controls
        
        ######Setting bounds for the lines not to exceed the number of slices when cutting###########
        self.horizontalLine.setBounds((0,self.ArrayDicom.shape[1]-1))
        self.verticalLine.setBounds((0,self.ArrayDicom.shape[0]-1))
        # self.obliqueLine.setBounds((0,self.ArrayDicom.shape[1]-1))
        self.horizontalLineSag.setBounds((0,self.sagittalView.shape[1]-1))
        self.verticalLineSag.setBounds((0,self.sagittalView.shape[0]-1))
        self.horizontalLineCor.setBounds((0,self.coronalView.shape[1]-1))
        self.verticalLineCor.setBounds((0,self.coronalView.shape[0]-1))

        ####################Setting Initial Values################
        self.horizontalLine.setValue(int(self.ArrayDicom.shape[1]/2)) # bring iso line above contrast controls
        self.verticalLine.setValue(int(self.ArrayDicom.shape[0]/2)) # bring iso line above contrast controls
        self.obliqueLine.setValue((int(self.ArrayDicom.shape[0]/2),int(self.ArrayDicom.shape[1]/2))) # bring iso line above contrast controls
        self.obliqueLineSag.setValue((int(self.sagittalView.shape[0]/2),int(self.sagittalView.shape[1]/2))) # bring iso line above contrast controls
        self.obliqueLineCor.setValue((int(self.coronalView.shape[0]/2),int(self.coronalView.shape[1]/2))) # bring iso line above contrast controls
        ####################When clicking get the view###############
        self.image.mousePressEvent = self.getPosAxial
        self.imageSagittal.mousePressEvent = self.getPosSagittal
        self.imageCoronal.mousePressEvent = self.getPosCoronal

        self.horizontalLineSag.setValue(int(self.sagittalView.shape[1]/2))
    

    ################################Function to generate the different views################
    def generateViews(self):
        # Step 1 3D numpy ---> VTK_ARRAY
        NumPy_data_shape = self.ArrayDicom.shape
        VTK_data = numpy_support.numpy_to_vtk(
                num_array=self.ArrayDicom.transpose(2,1,0).ravel(),  # ndarray contains the fitting result from the points. It is a 3D array
                deep=True,
                array_type=vtk.VTK_FLOAT)

        # Step 2 VTK_ARRAY ----> VTK__IMAGE_DATA
        img_vtk2 = vtk.vtkImageData()
        img_vtk2.GetPointData().SetScalars(VTK_data)
        img_vtk2.SetDimensions(NumPy_data_shape)
        img_vtk2.SetSpacing(self.ConstPixelSpacing[0], self.ConstPixelSpacing[1], self.ConstPixelSpacing[2])

        #######################origin of reslicing matrix################
        center= [0 + self.ConstPixelSpacing[0] * 0.5 * (0 + self.ArrayDicom.shape[0]),
                0 + self.ConstPixelSpacing[1] * 0.5 * (0 +  self.ArrayDicom.shape[1]),
                0 + self.ConstPixelSpacing[2] * 0.5 * (0 +  self.ArrayDicom.shape[2])]

        ####################coronal reslicing matrix####################
        coronal = vtk.vtkMatrix4x4()
        coronal.DeepCopy((1, 0, 0,0,
                        0, 0, 1, 0,
                        0,-1, 0, 1,
                        0, 0, 0, 1))
        ####################sagittal reslicing matrix####################
        sagittal = vtk.vtkMatrix4x4()
        sagittal.DeepCopy((0, 0,-1, center[0],
                        1, 0, 0, center[1],
                        0,-1, 0, center[2],
                        0, 0, 0, 1))
        ####################oblique basic reslicing matrix####################
        oblique = vtk.vtkMatrix4x4()
        oblique.DeepCopy((0.0, -0.515038, 0.857167, center[0],
                            0.0,  0.857167, 0.515038, center[1],
                            -1.0,  0.0,      0.0, center[2],
                            0.0,  0.0,      0.0,  1.0))

        self.sagittalView = self.viewOption(sagittal, img_vtk2)
        self.coronalView = self.viewOption(coronal, img_vtk2)
        self.obliqueView = self.ObliqueReconstruct(self.ArrayDicom, 45)
        self.sagittalView=np.flip(self.sagittalView, axis=0)
        self.sagittalView=self.sagittalView[:,:,::-1]

    def ObliqueReconstruct(self,Array,Angle,PlaneAngle=90):
        
        Array2=Array.copy()
        if(Angle!=0):
            Array2 = ndimage.rotate(Array2, angle=Angle, reshape=True)

        center= [0 + self.ConstPixelSpacing[0] * 0.5 * (0 + Array.shape[0]),
        0 + self.ConstPixelSpacing[1] * 0.5 * (0 +  Array.shape[1]),
        0 + self.ConstPixelSpacing[2] * 0.5 * (0 +  Array.shape[2])]


        # Step 1 3D numpy ---> VTK_ARRAY
        NumPy_data_shape = Array2.shape
        VTK_data = numpy_support.numpy_to_vtk(
                num_array=Array2.transpose(2,1,0).ravel(),  # ndarray contains the fitting result from the points. It is a 3D array
                deep=True,
                array_type=vtk.VTK_FLOAT)

        # Step 2 VTK_ARRAY ----> VTK__IMAGE_DATA
        img_vtk2 = vtk.vtkImageData()
        img_vtk2.GetPointData().SetScalars(VTK_data)
        img_vtk2.SetDimensions(NumPy_data_shape)
        img_vtk2.SetSpacing(self.ConstPixelSpacing[0], self.ConstPixelSpacing[1], self.ConstPixelSpacing[2])
        
        reslice = vtk.vtkImageReslice()
        reslice.SetInputData(img_vtk2)
        reslice.SetOutputDimensionality(3)

        # reslice.SetOutputSpacing(img_vtk2.GetSpacing())

        # *************************** #

        # reslice.SetResliceAxes(vtkResliceMatrix)

        
        reslice.SetResliceAxesOrigin(center[0], center[1],center[2])
		
        reslice.SetResliceAxesDirectionCosines(1, 0,  0, 
										0, mathh.cos(mathh.radians(PlaneAngle)),  -mathh.sin(mathh.radians(PlaneAngle)),
										0, mathh.sin(mathh.radians(PlaneAngle)),  mathh.cos(mathh.radians(PlaneAngle)) )

        reslice.SetInterpolationModeToLinear()


        reslice.Update() 

        reslicedImg = reslice.GetOutput() 


        reslicedNpImg = self.vtkToNumpy(reslicedImg)
        for i in range(2):
            reslicedNpImg=np.rot90(reslicedNpImg,1,axes=(0, 1))
        # reslicedNpImg = rotate(reslicedNpImg, angle=45)
        return reslicedNpImg



    #############applying the matrix on the volume to get the view##################
    def viewOption(self, vtkResliceMatrix, img_vtk2):
        # Extract a slice in the desired orientation
        reslice = vtk.vtkImageReslice()
        reslice.SetInputData(img_vtk2)
        reslice.SetOutputDimensionality(3)

        # reslice.SetOutputSpacing(img_vtk2.GetSpacing())

        # *************************** #

        reslice.SetResliceAxes(vtkResliceMatrix)
        reslice.SetInterpolationModeToLinear()


        reslice.Update() 

        reslicedImg = reslice.GetOutput() 


        reslicedNpImg = self.vtkToNumpy(reslicedImg) 
        
        return reslicedNpImg

    ################################Handles the scrolling of the mouse####################
    def eventFilter(self, watched, event):
        """
        Cases:
        Scrolling in axial : changing the red line in sagittal and coronal
        Scrolling in coronal: changing the red line in axial and blue line in sagittal 
        Scrolling in sagittal: changing the blue line in axial and red line in coronal
        """
        if (watched == self.ui.original.viewport() and 
            event.type() == QtCore.QEvent.Wheel):
            if event.angleDelta().y() > 0:
                self.indexOriginal+=1
                if(self.indexOriginal>self.ArrayDicom.shape[2]-1):
                    self.indexOriginal=self.ArrayDicom.shape[2]-1
                else:
                    img=self.ArrayDicom[:,:,self.indexOriginal]  
                    self.image = pg.ImageItem(img)      
                    self.ui.original.addItem(self.image)
                    self.horizontalLineSag.setValue((self.sagittalView.shape[1]-1)-self.indexOriginal)
                    # self.horizontalLineCor.setValue((self.coronalView.shape[1]-1)-self.indexOriginal)
            else:
                self.indexOriginal-=1
                if(self.indexOriginal<0):
                    self.indexOriginal=0
                else:
                    img=self.ArrayDicom[:,:,self.indexOriginal]  
                    self.image = pg.ImageItem(img)      
                    self.ui.original.addItem(self.image)
                    self.horizontalLineSag.setValue((self.ArrayDicom.shape[2])-self.indexOriginal)
                    self.horizontalLineCor.setValue((self.ArrayDicom.shape[2])-self.indexOriginal)

            return True
        elif (watched == self.ui.Sagittal.viewport() and 
            event.type() == QtCore.QEvent.Wheel):
            if event.angleDelta().y() > 0:
                self.indexSagittal+=1
                if(self.indexSagittal>self.sagittalView.shape[2]-1):
                    self.indexSagittal=self.sagittalView.shape[2]-1
                else:
                    img=self.sagittalView[:,:,self.indexSagittal]  
                    self.imageSagittal = pg.ImageItem(img)      
                    self.ui.Sagittal.addItem(self.imageSagittal)
                    self.verticalLine.setValue(self.indexSagittal)
            else:
                self.indexSagittal-=1
                if(self.indexSagittal<0):
                    self.indexSagittal=0
                else:
                    img=self.sagittalView[:,:,self.indexSagittal]  
                    self.imageSagittal = pg.ImageItem(img)      
                    self.ui.Sagittal.addItem(self.imageSagittal)
                    self.verticalLine.setValue(self.indexSagittal)
            return True
        elif (watched == self.ui.Coronal.viewport() and 
            event.type() == QtCore.QEvent.Wheel):
            if event.angleDelta().y() > 0:
                self.indexCoronal+=1
                if(self.indexCoronal>self.coronalView.shape[2]-1):
                    self.indexCoronal=self.coronalView.shape[2]-1
                else:
                    img=self.coronalView[:,:,self.indexCoronal]  
                    self.imageCoronal = pg.ImageItem(img)      
                    self.ui.Coronal.addItem(self.imageCoronal)
                    self.horizontalLine.setValue(self.indexCoronal)
            else:
                self.indexCoronal-=1
                if(self.indexCoronal<0):
                    self.indexCoronal=0
                else:
                    img=self.coronalView[:,:,self.indexCoronal]  
                    self.imageCoronal = pg.ImageItem(img)      
                    self.ui.Coronal.addItem(self.imageCoronal)
                    self.horizontalLine.setValue(self.indexCoronal)
            return True
        elif (watched == self.ui.Oblique.viewport() and 
            event.type() == QtCore.QEvent.Wheel):
            if event.angleDelta().y() > 0:
                self.indexOblique+=1
                img=self.obliqueView[:,:,self.indexOblique]  
                self.image = pg.ImageItem(img)      
                self.ui.Oblique.addItem(self.image)
                if self.whichview==1:
                    self.obliqueLine.setValue((self.indexOblique*(512/self.obliqueView.shape[2]),self.indexOblique*(512/self.obliqueView.shape[2]))) # bring iso line above contrast controls
                elif self.whichview==2:
                    self.obliqueLineSag.setValue((self.obliqueView.shape[2]-(self.indexOblique*(self.ArrayDicom.shape[0]/self.obliqueView.shape[2])),self.obliqueView.shape[2]-(self.indexOblique*(self.ArrayDicom.shape[2]/self.obliqueView.shape[2])))) # bring iso line above contrast controls
                elif self.whichview==3:
                    self.obliqueLineCor.setValue((self.obliqueView.shape[2]-(self.indexOblique*(self.ArrayDicom.shape[0]/self.obliqueView.shape[2])),self.obliqueView.shape[2]-(self.indexOblique*(self.ArrayDicom.shape[2]/self.obliqueView.shape[2])))) # bring iso line above contrast controls

            else:
                self.indexOblique-=1
                if(self.indexOblique<0):
                    self.indexOblique=0
                else:
                    img=self.obliqueView[:,:,self.indexOblique]  
                    self.image = pg.ImageItem(img)      
                    self.ui.Oblique.addItem(self.image)
                    if self.whichview==1:
                        self.obliqueLine.setValue((self.indexOblique*(512/self.obliqueView.shape[2]),self.indexOblique*(512/self.obliqueView.shape[2]))) # bring iso line above contrast controls
                    elif self.whichview==2:
                        self.obliqueLineSag.setValue((self.obliqueView.shape[2]-(self.indexOblique*(self.ArrayDicom.shape[0]/self.obliqueView.shape[2])),self.obliqueView.shape[2]-(self.indexOblique*(self.ArrayDicom.shape[2]/self.obliqueView.shape[2])))) # bring iso line above contrast controls
                    elif self.whichview==3:
                        self.obliqueLineCor.setValue((self.obliqueView.shape[2]-(self.indexOblique*(self.ArrayDicom.shape[0]/self.obliqueView.shape[2])),self.obliqueView.shape[2]-(self.indexOblique*(self.ArrayDicom.shape[2]/self.obliqueView.shape[2])))) # bring iso line above contrast controls

 
            return True

        return super().eventFilter(watched, event)

def main():

    app = QtWidgets.QApplication(sys.argv)
    application = ApplicationWindow()
    application.show()
    app.exec_()
    
    

if __name__ == '__main__':
    main()

