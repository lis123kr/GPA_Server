from django.shortcuts import render
from django.conf import settings
import pandas as pd
from . import Excel
# Create your views here.
xlsx = object()
filename = None
def index(request):
	import os
	for f, _, files in os.walk(settings.MEDIA_ROOT):
		for file in files:
			os.remove(f+'\\'+file)
	return render(request, 'blog/index.html', {})

def gpa_xlsx(request):
	global xlsx, filename
	filename = request.FILES.get('uploaded_file')
	if filename is not None:
		# start_time = time.time()
		xlsx = pd.ExcelFile(filename)
		sheet_names = xlsx.sheet_names
		sheets = list()
		for x in sheet_names:
			sheets.append(xlsx.parse(x))
		cols = sheets[0].columns

		return render(request, 'blog/gpa_xlsx.html', {
			'sheet_names': sheet_names,
			'filename' : str(filename),
			'cols': cols,
		})
	else:
		return render(request, 'blog/gpa_xlsx.html', {})
def replace_(L):
	if L is not None:
	    L = L.replace('\r\n',',').replace('\r',',').replace('\n',',').replace(',,',',').replace(',,',',').replace(',,',',').replace(',,',',').replace(',,',',').strip(',').split(',')
	    if L[-1] == '':
	        L.pop()
	    if L[0] == '':
	        L.pop(0)
	    return L

def some_streaming_xlsx_view(request):
	global xlsx, filename
	selected_sheets = request.POST.getlist('_selected_')
	t_genome = replace_(request.POST.get('genome'))
	t_region = replace_(request.POST.get('region'))
	t_orf = replace_(request.POST.get('orf'))
	t_ncr = replace_(request.POST.get('ncr'))

	c_dm_position = request.POST.get('dm_position')
	c_dm_seq = request.POST.get('dm_seq')
	c_seq = request.POST.get('seq')
	c_genome = str(request.POST.get('col_genome'))
	c_region = str(request.POST.get('col_region'))
	c_orf = request.POST.get('col_orf')

	atype = request.POST.get('atype')

	excel = Excel.Excel(xlsx, str(filename), selected_sheets, 
                    c_dm_position, c_dm_seq, c_genome, c_region, c_orf, c_seq, 
                    t_genome, t_region, t_orf, t_ncr)
	import os, zipfile, io
	from django.http import HttpResponse
	path=os.path.join(settings.MEDIA_ROOT, str(filename))
	p = excel.Analyze(settings.MEDIA_ROOT, atype)

	s = io.BytesIO()
	zf = zipfile.ZipFile(s, 'w')
	for folder, _, files in os.walk(settings.MEDIA_ROOT):
		for file in files:
			print(os.path.relpath(os.path.join(folder,file), settings.MEDIA_ROOT))
			zf.write(os.path.join(folder, file), os.path.relpath(os.path.join(folder,file), settings.MEDIA_ROOT), compress_type = zipfile.ZIP_DEFLATED)
	zf.close()

	response = HttpResponse(s.getvalue(), content_type="application/x-zip-compressed")
	response['Content-Disposition'] = 'attachment; filename='+str(filename)+'.zip'
	return response

