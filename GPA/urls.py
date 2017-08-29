from django.conf.urls import url, include
from django.contrib import admin
from . import views

urlpatterns = [
    url(r'^$', views.index, name='index'),
    url(r'gpa_xlsx$', views.gpa_xlsx, name='gpa_xlsx'),
    url(r'save$', views.some_streaming_xlsx_view, name='some_view'),
]
