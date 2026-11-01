import os
from chimerax.core.tools import ToolInstance
from chimerax.ui import MainToolWindow
from Qt.QtWidgets import (QWidget,QVBoxLayout,QFormLayout,QComboBox,QDoubleSpinBox,QSpinBox,
                          QPushButton,QLabel,QLineEdit,QHBoxLayout,QGroupBox,QCheckBox)
from .engine import run_docking

class CholesterolDockingTool(ToolInstance):
    SESSION_ENDURING=False
    SESSION_SAVE=False
    def __init__(self,session,tool_name):
        super().__init__(session,tool_name)
        self.tool_window=MainToolWindow(self)
        root=QWidget(); layout=QVBoxLayout(root)
        title=QLabel('<b>Automated Cholesterol Docking V7</b><br>V6.1 membrane detector + neighbor optimization + confidence diagnostics')
        title.setWordWrap(True);layout.addWidget(title)
        form=QFormLayout();layout.addLayout(form)
        self.map_combo=QComboBox();self.protein_combo=QComboBox();form.addRow('Cryo-EM map',self.map_combo);form.addRow('Protein model',self.protein_combo)
        self.detect=QDoubleSpinBox();self.detect.setDecimals(4);self.detect.setRange(-10,100);self.detect.setValue(.141);self.detect.setSingleStep(.01);form.addRow('Detection contour',self.detect)
        self.support=QDoubleSpinBox();self.support.setDecimals(4);self.support.setRange(-10,100);self.support.setValue(.05);self.support.setSingleStep(.01);form.addRow('Support level',self.support)
        self.quality=QComboBox();self.quality.addItems(['Fast (500)','Standard (1500)','Thorough (5000)']);self.quality.setCurrentIndex(1);form.addRow('Search quality',self.quality)
        self.maxsites=QSpinBox();self.maxsites.setRange(1,500);self.maxsites.setValue(60);form.addRow('Maximum sites',self.maxsites)
        self.out=QLineEdit(os.path.expanduser('~/cholesterol_v7'));form.addRow('Output folder',self.out)
        adv=QGroupBox('Advanced');af=QFormLayout(adv);layout.addWidget(adv)
        self.res=QDoubleSpinBox();self.res.setRange(3,20);self.res.setValue(6);af.addRow('Protein map resolution (A)',self.res)
        self.smooth=QDoubleSpinBox();self.smooth.setRange(.5,10);self.smooth.setValue(3);af.addRow('Blob smoothing sigma (A)',self.smooth)
        self.blob=QDoubleSpinBox();self.blob.setDecimals(3);self.blob.setRange(.01,.9);self.blob.setValue(.10);af.addRow('Blob threshold fraction',self.blob)
        self.pmask=QDoubleSpinBox();self.pmask.setRange(0,6);self.pmask.setValue(2.5);af.addRow('Residual protein mask (A)',self.pmask)
        row=QHBoxLayout();self.refresh=QPushButton('Refresh models');self.run_btn=QPushButton('Find & Dock Cholesterols');row.addWidget(self.refresh);row.addWidget(self.run_btn);layout.addLayout(row)
        self.status=QLabel('Ready. Results include HIGH/MEDIUM/LOW confidence models and site_diagnostics.csv.');self.status.setWordWrap(True);layout.addWidget(self.status)
        self.refresh.clicked.connect(self._refresh);self.run_btn.clicked.connect(self._run);self._refresh()
        self.tool_window.ui_area.setLayout(QVBoxLayout());self.tool_window.ui_area.layout().addWidget(root);self.tool_window.manage('side')

    def _refresh(self):
        self.map_combo.clear();self.protein_combo.clear()
        for m in self.session.models.list():
            if hasattr(m,'data') and hasattr(m,'interpolated_values'):
                self.map_combo.addItem('#%s  %s'%(m.id_string,m.name),m)
            try:
                if hasattr(m,'atoms') and hasattr(m,'residues') and len(m.atoms)>0:
                    self.protein_combo.addItem('#%s  %s'%(m.id_string,m.name),m)
            except Exception:pass

    def _run(self):
        m=self.map_combo.currentData();p=self.protein_combo.currentData()
        if m is None or p is None:self.status.setText('Select one map and one protein model.');return
        orientations=[500,1500,5000][self.quality.currentIndex()]
        out=os.path.abspath(os.path.expanduser(self.out.text()))
        self.status.setText('Running docking. Progress is reported in the ChimeraX Log...')
        try:
            result=run_docking(self.session,[m],[p],out,detection_level=self.detect.value(),support_level=self.support.value(),orientations=orientations,maxsites=self.maxsites.value(),protein_map_resolution=self.res.value(),blob_smoothing=self.smooth.value(),blob_threshold_fraction=self.blob.value(),protein_mask_radius=self.pmask.value())
            a=result['accepted'];h=sum(x['confidence']=='HIGH' for x in a);med=sum(x['confidence']=='MEDIUM' for x in a);lo=sum(x['confidence']=='LOW' for x in a)
            self.status.setText('Complete: %d accepted — HIGH %d, MEDIUM %d, LOW %d. Diagnostics: %s'%(len(a),h,med,lo,os.path.join(out,'site_diagnostics.csv')))
        except Exception as e:
            self.status.setText('Failed: %s (see ChimeraX Log)'%e);raise
