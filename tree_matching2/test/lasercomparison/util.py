# -*- coding: utf-8 -*-
from openalea.mtg.io import *
import cPickle as pickle


def writefile(fn, obj):
    f = open(fn,'wb')
    pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)
    f.close()
    
def readfile(fn, mode='rb'):
    f = open(fn,mode)
    obj = pickle.load(f)
    f.close()
    return obj

def checkfile(fn):
    import os.path
    return os.path.exists(fn)
    
    
def writeMTGfile(fn, g, properties=[('XX','REAL'), ('YY','REAL'), ('ZZ','REAL'), ('radius','REAL')], nb_tab=20):
    
    if properties == []:
      properties = [(p, 'REAL') for p in g.property_names() if p not in ['edge_type', 'index', 'label']]
    str = write_mtg(g, properties, nb_tab=nb_tab)
    f = open(fn, 'w')
    f.write(str)
    f.close()
  

def convertToStdMTG(g):
  from copy import deepcopy
  newg = deepcopy(g)
  
  pdic = newg.property('position')
  xx = {}
  yy = {}
  zz = {}
  r = {}
  
  for i,v in pdic.iteritems():
      xx[i]=v.x
      yy[i]=v.y
      zz[i]=v.z
      r[i]=0
  
  newg.add_property('XX')
  newg.add_property('YY')
  newg.add_property('ZZ')
  newg.add_property('radius')
  newg.property('XX').update(xx)
  newg.property('YY').update(yy)
  newg.property('ZZ').update(zz)
  newg.property('radius').update(r)
  del newg.properties()['position']
  return newg
  """
  from openalea.mtg import MTG
  
  
  
  def addProperties(mtg, vid, position, radius):
      mtg.property('XX')[vid] = position.x
      mtg.property('YY')[vid] = position.y
      mtg.property('ZZ')[vid] = position.z
      mtg.property('radius')[vid] = radius
      
  mtg = MTG()
  mtg.add_property('XX')
  mtg.add_property('YY')
  mtg.add_property('ZZ')
  mtg.add_property('radius')
    
  plantroot = mtg.root
  branchroot = mtg.add_component(plantroot,label='B')
  noderoot = mtg.add_component(branchroot,label='N')
  
  pdic = g.property('position')
  
  for k,v in pdic.iteritems():
    parentid = g.parent(k)
    r = g.property('radius')[k]
    
    if parentid == None:
      addProperties(mtg, k, v, r)
    else:
      label = g.label(k)
      if label == 'N':
        vid = mtg.add_child(parentid,edge_type='<',label='N')
      else:
        vid = mtg.add_child(parentid,edge_type='+',label='B')
        
      addProperties(mtg, vid, v, r)
    
  return mtg
  """

def convertToMyMTG(g):
  from openalea.mtg import MTG
  
  def addProperties(mtg, vid, px, py, pz, radius):
      mtg.property('position')[vid] = Vector3(px,py,pz)
      mtg.property('radius')[vid] = radius
      
  mtg = MTG()
  mtg.add_property('position')
  mtg.add_property('radius')
    
  plantroot = mtg.root
  branchroot = mtg.add_component(plantroot,label='B')
  noderoot = mtg.add_component(branchroot,label='N')
  
  rdic = g.property('radius')
  
  for k,r in rdic.iteritems():
    parentid = g.parent(k)
    px = g.property('XX')[k]
    py = g.property('YY')[k]
    pz = g.property('ZZ')[k]
    
    if parentid == None:
      addProperties(mtg, k, px, py, pz, r)
    else:
      label = g.label(k)
      if label == 'N':
        vid = mtg.add_child(parentid,edge_type='<',label='N')
      else:
        vid = mtg.add_child(parentid,edge_type='+',label='B')
        
      addProperties(mtg, vid, px, py, pz, r)
    
  return mtg
  
  
  