from django.urls import path

from . import views

urlpatterns = [
    path('', views.index, name='index'),
    path('regular', views.maxcut_sim, name='maxcut_sim'),
    path('get-maxcut', views.get_maxcut, name='get_maxcut'),
    path('get-mis', views.get_mis, name='get_mis'),
    path('init-mc-jobs', views.init_mc_jobs, name='init_mc_jobs'),
    path('init-mis-jobs', views.init_mis_jobs, name='init_mis_jobs'),
    path('companyfleets', views.mis_sim, name='mis_sim'),
]
