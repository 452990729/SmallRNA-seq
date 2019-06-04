
$(function() {
	var b,c, d, e,f;
	if (document.getElementById("eta")){
		e=echarts.init(document.getElementById("etll"));
		e.setOption({
		legend: {
			data: names
		},
		grid:[{height:0,top:'100%'}],
		toolbox: {
			feature:{
				myTool1:{
					show:1,
					title:'显示全部',
					icon:'path://M300,120L107.143,300L0,200V80l107.143,99.688L300,0V120z',
					onclick:function(){
legendallopen(e);
legendallopen(ate);
legendallopen(bte);
					}
				},myTool2:{
					show:1,
					title:'取消全部',
					icon:'path://M300,300h-42.894L150,194.196L44.221,300H0v-42.857L105.779,150L0,44.196V0h44.221L150,105.804L257.106,0H300v44.196L194.221,150L300,255.804V300z',
		onclick:function(){
legendallclose(e);
legendallclose(ate);
legendallclose(bte);
					}
				}
			}
			
		},
		xAxis: [{
			type: "category",show:0,
			data: ["1"]
		}],
		yAxis: [{
			type: "value",show:0
		}],
				color:['#058DC7','#50B432','#e39475','#DDDF00','#24CBE5','#64E572','#FF9655','#6AF9C4','#AF94FF','#FFF263','#05C766','#10B6Ff','#43a102','#065fb9','#a2b700','#ffc45c','#b4f593','#5aacdc','#4acca0','#FF6633','#61a2ff','#f17D54'],
				series: legendata})};
		if (document.getElementById("etu")){
			f=echarts.init(document.getElementById("etlu"));
	f.setOption({
		legend: {
			data: names
		},
		grid:[{height:0,top:'100%'}],
		toolbox: {
			feature:{
				myTool1:{
					show:1,
					title:'显示全部',
					icon:'path://M300,120L107.143,300L0,200V80l107.143,99.688L300,0V120z',
					onclick:function(){
legendallopen(f);
legendallopen(ute);
					}
				},myTool2:{
					show:1,
					title:'取消全部',
					icon:'path://M300,300h-42.894L150,194.196L44.221,300H0v-42.857L105.779,150L0,44.196V0h44.221L150,105.804L257.106,0H300v44.196L194.221,150L300,255.804V300z',
		onclick:function(){
legendallclose(f);
legendallclose(ute);
					}
				}
			}
			
		},
		xAxis: [{
			type: "category",show:0,
			data: ["1"]
		}],
		yAxis: [{
			type: "value",show:0
		}],
		color:['#058DC7','#50B432','#e39475','#DDDF00','#24CBE5','#64E572','#FF9655','#6AF9C4','#AF94FF','#FFF263','#05C766','#10B6Ff','#43a102','#065fb9','#a2b700','#ffc45c','#b4f593','#5aacdc','#4acca0','#FF6633','#61a2ff','#f17D54'],
		series: legendata})};
	if(document.getElementById("etu")){
		ute = echarts.init(document.getElementById("etu"));
	ute.setOption({
		legend: {
			data: names,
			show:0
		},title:{show:0,text:"每条染色体的覆盖深度和覆盖率"},
		tooltip:{show:1,trigger:'item',formatter: "{a}</br>{b} : {c}"},
		toolbox: {
			show: 1,
			itemSize:30,
			feature: {
				restore:{show:1},
				saveAsImage: {
					show: 1
				}
				
			},
			iconStyle:{emphasis:{borderColor:'#0f5b7c'}}
		},
		dataZoom:[{show:1}],
		xAxis: [{
			type: "category",
			data: ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"]
		}],
		yAxis: [{
			type: "value",
			name: "覆盖深度",
			splitArea: {
				show: !0
			},
			axisLabel: {
				formatter: "{value}x"
			},
			min: 0
		}, {
			type: "value",
			name: "覆盖率",
			splitArea: {
				show: !0
			},
			axisLabel: {
				formatter: function(value){var data=value*100;return data.toString()+"%"}
			},
			min: 0,
			max: 1
		}],
		color:['#058DC7','#058DC7','#50B432','#50B432','#e39475','#e39475','#DDDF00','#DDDF00','#24CBE5','#24CBE5','#64E572','#64E572','#FF9655','#FF9655','#6AF9C4','#6AF9C4','#AF94FF','#AF94FF','#FFF263','#FFF263','#05C766','#05C766','#10B6Ff','#10B6Ff','#43a102','#43a102','#065fb9','#065fb9','#a2b700','#a2b700','#ffc45c','#ffc45c','#b4f593','#b4f593','#5aacdc','#5aacdc','#4acca0','#4acca0','#FF6633','#FF6633','#61a2ff','#61a2ff','#f17D54','#f17D54'],
		series: chrdata
})}; 
if(document.getElementById("eta")){
	ate = echarts.init(document.getElementById("eta")), 
	ate.setOption({
		legend: {
			data: names,
show:0
		},title:{show:0,text:"不同测序深度的碱基比例"},
		toolbox: {
			show: 1,
			itemSize:30,
			feature: {
				saveAsImage: {
					show: 1
				}
			},
			iconStyle:{emphasis:{borderColor:'#0f5b7c'}}
		},tooltip:{show:1,trigger:'item',formatter: "{a}"},
		xAxis: [{
			type: "value",
			name: "深度分布",
			max: depthavg
		}],
		yAxis: [{
			type: "value",
			splitArea: {
				show: !0
			},
			min: 0,
			axisLabel: {
				formatter: function(value){var data=(value*100).toFixed(2);return data.toString()+"%"}
			}
		}],
		color:['#058DC7','#50B432','#e39475','#DDDF00','#24CBE5','#64E572','#FF9655','#6AF9C4','#AF94FF','#FFF263','#05C766','#10B6Ff','#43a102','#065fb9','#a2b700','#ffc45c','#b4f593','#5aacdc','#4acca0','#FF6633','#61a2ff','#f17D54'],
		series: depthdist
	})}; 
	if(document.getElementById("etb")){
	bte = echarts.init(document.getElementById("etb")), 
	bte.setOption({
		legend: {
			data: names,
show:0
		},title:{show:0,text:"不同深度上的累积碱基比例"},
		toolbox: {
			show: 1,
			itemSize:30,
			feature: {
				saveAsImage: {
					show: 1
				}
			},
			iconStyle:{emphasis:{borderColor:'#0f5b7c'}}
		},tooltip:{show:1,trigger:'item',formatter: "{a}"},
		xAxis: [{
			type: "value",
			name: "累积深度",
			max: depthavg
		}],
		yAxis: [{
			type: "value",
			splitArea: {
				show: !0
			},
			axisLabel: {
				formatter: function(value){var data=value*100;return data.toString()+"%"}
			},
			min: 0,
			max: 1
		}],
		color:['#058DC7','#50B432','#e39475','#DDDF00','#24CBE5','#64E572','#FF9655','#6AF9C4','#AF94FF','#FFF263','#05C766','#10B6Ff','#43a102','#065fb9','#a2b700','#ffc45c','#b4f593','#5aacdc','#4acca0','#FF6633','#61a2ff','#f17D54'],
		series: depthcumu
	})}; 
	if(document.getElementById("etu")){echarts.connect([ute,f])};
	if(document.getElementById("eta")){echarts.connect([ate,e])};
	if(document.getElementById("etb")){echarts.connect([bte,e])};
	if(document.getElementById("QMh")){
	qmh=echarts.init(document.getElementById("QMh")), 
	qmh.setOption({
		title: {
			text: "Classification of Raw Reads"
		},
		toolbox: {
			show: 1,
			itemSize:30,
			feature: {
				saveAsImage: {
					show: 1
				}
			},
			iconStyle:{emphasis:{borderColor:'#0f5b7c'}}
		},
		tooltip : {
        trigger: 'item',
        formatter: "{b} : {c} ({d}%)"
    },series:[{type:'pie',label:{normal:{textStyle:{color:'#000000'},formatter: "{b} : {c} ({d}%)"}},radius:'60%',data:nbdata2}]
	})
	;a = document.getElementById("qmname"), 
	a.onchange = function() {
		var c,b = a.options[document.getElementById("qmname").options.selectedIndex].value;
		pqm=document.getElementById("pqm"),
		nqm=document.getElementById("nqm"),
		nbdata2 = datanow2[b], qmh.setOption({series:[{type:'pie',data:nbdata2}]});
	}, $("#pqm").click(function() {
		a.selectedIndex = a.selectedIndex - 1, nqm.style.visibility="visible",-1 == a.selectedIndex && (a.selectedIndex = $("#qmname").find("option:first").index(),pqm.style.visibility="hidden"), a.onchange();
	}), $("#nqm").click(function() {
		a.selectedIndex = a.selectedIndex + 1,  pqm.style.visibility="visible",-1 == a.selectedIndex && (a.selectedIndex = $("#qmname").find("option:last").index(),nqm.style.visibility="hidden"), a.onchange();
})};  
if(document.getElementById("ERm")){
	erm=echarts.init(document.getElementById("ERm")),
	erm.setOption({
		title: {
			text: "Error rate distribution along reads"
		},
		tooltip:{show:1,trigger:'axis',formatter: "{b} : {c}"},
		toolbox: {
			show: 1,
			itemSize:30,
			feature: {
				saveAsImage: {
					show: 1
				}
			},
			iconStyle:{emphasis:{borderColor:'#0f5b7c'}}
		},
		xAxis:[
		{
			type:'category',
			data:num,
			axisLabel:{interval:49}
		}
		],yAxis:[{type:'value',axisLabel: {formatter: "{value}%"}}],series:[{type:'line',areaStyle:{normal:{}},data:ermdata}]
}),b = document.getElementById("ermname"), 
	b.onchange = function() {
		ermnow = b.options[document.getElementById("ermname").options.selectedIndex].value, 
		perm=document.getElementById("perm"),
		nerm=document.getElementById("nerm"),
		ermdata = ermd[ermnow], erm.setOption({series:[{data:ermdata}]})
	};$("#perm").click(function() {
		b.selectedIndex = b.selectedIndex - 1, nerm.style.visibility="visible",-1 == b.selectedIndex && (b.selectedIndex = $("#ermname").find("option:first").index(),perm.style.visibility="hidden"), b.onchange();
	}), $("#nerm").click(function() {
		b.selectedIndex = b.selectedIndex + 1, perm.style.visibility="visible",-1 == b.selectedIndex && (b.selectedIndex = $("#ermname").find("option:last").index(),nerm.style.visibility="hidden"), b.onchange();
	})};
if(document.getElementById("GCz")){
	gcz=echarts.init(document.getElementById("GCz")),
    gcz.setOption({
		title: {
			text: "Bases content along reads"
		},
		tooltip:{show:1,trigger:'axis',formatter: "{b0}<br />{a0}: {c0}<br />{a1}: {c1}<br />{a2}: {c2}<br />{a3}: {c3}<br />{a4}: {c4}"},
		toolbox: {
			show: 1,
			itemSize:30,
			feature: {
				saveAsImage: {
					show: 1
				}
			},
			iconStyle:{emphasis:{borderColor:'#0f5b7c'}}
		},xAxis:[{
			type:'category',
			data:num,
			axisLabel:{interval:49}
		}
		],yAxis:[{type:'value',min:0,max:50}],
		legend:{data:['A','T','G','C','N'],top:'bottom'},series:gczdata1
		})
	;c = document.getElementById("gczname"), 
	c.onchange = function() {
		gcznow = c.options[document.getElementById("gczname").options.selectedIndex].value, 
		pgcz=document.getElementById("pgcz"),
		ngcz=document.getElementById("ngcz"),
		gczdata1 = gczd1[gcznow], 
		gcz.setOption({series:gczdata1})
	};

	 $("#pgcz").click(function() {
		c.selectedIndex = c.selectedIndex - 1, ngcz.style.visibility="visible",-1 == c.selectedIndex && (c.selectedIndex = $("#gczname").find("option:first").index(),pgcz.style.visibility="hidden"), c.onchange();
	}), $("#ngcz").click(function() {
		c.selectedIndex = c.selectedIndex + 1, pgcz.style.visibility="visible",-1 == c.selectedIndex && (c.selectedIndex = $("#gczname").find("option:last").index(),ngcz.style.visibility="hidden"), c.onchange();
	})}; 
	if(document.getElementById("ERz")){
		erz=echarts.init(document.getElementById("ERz")),
		erz.setOption({
		title: {
			text: "Quality score distribution along reads"
		},
		tooltip:{show:1,trigger:'axis',formatter: "{b} : {c}"},
		toolbox: {
			show: 1,
			itemSize:30,
			feature: {
				saveAsImage: {
					show: 1
				}
			},
			iconStyle:{emphasis:{borderColor:'#0f5b7c'}}
		},
		xAxis:[
		{
			type:'category',
			data:num,
			axisLabel:{interval:49}
		}
		],yAxis:[{type:'value',min:0,max:50}],series:[{type:'line',data:erzdata}]
		})
		
	;d = document.getElementById("erzname"),
	d.onchange = function() {
		erznow = d.options[document.getElementById("erzname").options.selectedIndex].value, 
		perz=document.getElementById("perz"),
		nerz=document.getElementById("nerz"),
		erzdata = erzd[erznow], erz.setOption({
			series:[{data: erzdata}]
		});
	};$("#perz").click(function() {
		d.selectedIndex = d.selectedIndex - 1, nerz.style.visibility="visible", -1 == d.selectedIndex && (d.selectedIndex = $("#erzname").find("option:first").index(), perz.style.visibility="hidden"), d.onchange();
	}), $("#nerz").click(function() {
		d.selectedIndex = d.selectedIndex + 1, perz.style.visibility="visible", -1 == d.selectedIndex && (d.selectedIndex = $("#erzname").find("option:last").index(), nerz.style.visibility="hidden"), d.onchange();
	})};
	if(document.getElementById("SNv_exoni1")){
	snv_exoni1=echarts.init(document.getElementById("SNv_exoni1"))
	snv_exoni1.setOption({
		title: [{
			text: "编码区上不同类型的SNV数目分布",
			left:'left'
		},{
			text: "基因组不同区域上SNV数目分布",
			left:'50%'
		}],
		toolbox: {
			show: 1,
			itemSize:30,
			feature: {
				saveAsImage: {
					show: 1
				}
			},
			iconStyle:{emphasis:{borderColor:'#0f5b7c'}}
		},
		tooltip : {
        trigger: 'item',
        formatter: "{b} : {c} ({d}%)"
    },series:[
	{type:'pie',color:['#058DC7', '#50B432', '#e39475', '#DDDF00', '#24CBE5', '#64E572', '#FF9655', '#FFF263', '#6AF9C4'],
label:{normal:{textStyle:{color:'#000000'},formatter: "{b} : {c} ({d}%)"}},radius:'50%',center:['30%','50%'],data:snvexoni2},
	{type:'pie',color:['#058DC7', '#50B432', '#e39475', '#DDDF00', '#24CBE5', '#64E572', '#FF9655', '#FFF263', '#6AF9C4'],
label:{normal:{textStyle:{color:'#000000'},formatter: "{b} : {c} ({d}%)"}},radius:'50%',center:['70%','50%'],data:snvfunc2}]
	});
	snvpieb = document.getElementById("snvpie"), 

psnvp=document.getElementById("psnvp"),
nsnvp=document.getElementById("nsnvp"),
	snvpieb.onchange = function() {
		name = snvpieb.options[document.getElementById("snvpie").options.selectedIndex].value, 
		snvfunc2 = snvfuncz2[name], 
		snvexoni2 = snvexoniz2[name], 
		snv_exoni1.setOption({series:[{data:snvexoni2},{data:snvfunc2}]})
	};
	$("#psnvp").click(function() {
		snvpieb.selectedIndex = snvpieb.selectedIndex - 1,nsnvp.style.visibility="visible", -1 == snvpieb.selectedIndex && (snvpieb.selectedIndex = $("#snvpie").find("option:first").index(),psnvp.style.visibility="hidden"), snvpieb.onchange();
	}), $("#nsnvp").click(function() {
		snvpieb.selectedIndex = snvpieb.selectedIndex + 1,psnvp.style.visibility="visible", -1 == snvpieb.selectedIndex && (snvpieb.selectedIndex = $("#snvpie").find("option:last").index(),nsnvp.style.visibility="hidden"), snvpieb.onchange();
	})}; 
	if(document.getElementById("SNv_type")){
	snv_type=echarts.init(document.getElementById("SNv_type")),
	snv_type.setOption({
		title: {
			text: "基因组SNV特征"
		},
		toolbox: {
			show: 1,
			itemSize:30,
			feature: {
				saveAsImage: {
					show: 1
				}
			},
			iconStyle:{emphasis:{borderColor:'#0f5b7c'}}
		},
		tooltip : {
        trigger: 'item',
        formatter: "{b} : {c} ({d}%)"
    },series:[
	{type:'pie',
	label:{normal:{textStyle:{color:'#000000'},formatter: "{b} : {c} ({d}%)"}},
	radius:'35%',
	center : ['30%', '25%'],
	data:snp_homhet},
	{type:'pie',
	label:{normal:{textStyle:{color:'#000000'},formatter: "{b} : {c} ({d}%)"}},
	radius:'35%',
	center : ['70%', '25%'],
	data:snp_newold},
	{type:'pie',
	label:{normal:{textStyle:{color:'#000000'},formatter: "{b} : {c} ({d}%)"}},
	radius:'35%',
	center : ['30%', '75%'],
	data:snp_tstv},
	{type:'pie',
	label:{normal:{textStyle:{color:'#000000'},formatter: "{b} : {c} ({d}%)"}},
	radius:'35%',
	center : ['70%', '75%'],
	data:snp_newtstv}
	]
});
snvhunb = document.getElementById("snvhun"), 
psnvhun=document.getElementById("psnvhun"),
nsnvhun=document.getElementById("nsnvhun"),
snvhunb.onchange = function() {
		name = snvhunb.options[document.getElementById("snvhun").options.selectedIndex].value, 
		snp_homhet = snp_homhetz[name], 
		snp_newold = snp_newoldz[name], 
		snp_tstv=snp_tstvz[name],
		snp_newtstv=snp_newtstvz[name],
		snv_type.setOption({
			series:[
	{data:snp_homhet},	{data:snp_newold},	{data:snp_tstv},	{data:snp_newtstv}	]
	})
	};
	$("#psnvhun").click(function() {
		snvhunb.selectedIndex = snvhunb.selectedIndex - 1,nsnvhun.style.visibility="visible", -1 == snvhunb.selectedIndex && (snvhunb.selectedIndex = $("#snvhun").find("option:first").index(),psnvhun.style.visibility="hidden"), snvhunb.onchange();
	}), $("#nsnvhun").click(function() {
		snvhunb.selectedIndex = snvhunb.selectedIndex + 1,psnvhun.style.visibility="visible", -1 == snvhunb.selectedIndex && (snvhunb.selectedIndex = $("#snvhun").find("option:last").index(),nsnvhun.style.visibility="hidden"), snvhunb.onchange();
	})};
	if(document.getElementById("INdellen")){
	indellen=echarts.init(document.getElementById("INdellen"))
	indellen.setOption({
		title: {
			text: "InDel长度分布"
		},
		legend:{data:['cds','ncds'],top:'bottom'},
		tooltip:{show:1,trigger:'axis',formatter: "{b} : {c}"},
		toolbox: {
			show: 1,
			itemSize:30,
			feature: {
				saveAsImage: {
					show: 1
				}
			},
			iconStyle:{emphasis:{borderColor:'#0f5b7c'}}
		},
		xAxis:[
		{
			type:'category',
			data:[-30, -29, -28, -27, -26, -25, -24, -23, -22, -21, -20, -19, -18, -17, -16, -15, -14, -13, -12, -11, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30],
			axisLabel:{interval:2}
		}
		],yAxis:[{type:'value'}],
		series: [{
			type:'bar',
			name: 'cds',
			data: indel_len
		}, {
			type:'bar',
			name: 'ncds',
			data: indel_lenn
		}]
});
indellb = document.getElementById("indell"), 
pidl=document.getElementById("pidl"),
nidl=document.getElementById("nidl"),
indellb.onchange = function() {
		name = indellb.options[document.getElementById("indell").options.selectedIndex].value, 
		indel_len = indel_lenz[name], 
		indel_lenn = indel_lennz[name], 
		indellen.setOption({series:[{data:indel_len},{data:indel_lenn}]})
	};
	$("#pidl").click(function() {
		indellb.selectedIndex = indellb.selectedIndex - 1,nidl.style.visibility="visible", -1 == indellb.selectedIndex && (indellb.selectedIndex = $("#indell").find("option:first").index(),pidl.style.visibility="hidden"), indellb.onchange();
	}), $("#nidl").click(function() {
		indellb.selectedIndex = indellb.selectedIndex + 1,pidl.style.visibility="visible", -1 == indellb.selectedIndex && (indellb.selectedIndex = $("#indell").find("option:last").index(),nidl.style.visibility="hidden"), indellb.onchange();
	})};
	
	if(document.getElementById("INdel_exoni1")){
	indel_exoni1=echarts.init(document.getElementById("INdel_exoni1"))
	indel_exoni1.setOption({
		title: [{
			text: "编码区上不同类型的INDEL数目分布",
			left:'left'
		},{
			text: "基因组不同区域上InDel数目分布",
			left:'50%'
		}],
		toolbox: {
			show: 1,
			itemSize:30,
			feature: {
				saveAsImage: {
					show: 1
				}
			},
			iconStyle:{emphasis:{borderColor:'#0f5b7c'}}
		},
		tooltip : {
        trigger: 'item',
        formatter: "{b} : {c} ({d}%)"
    },series:[{type:'pie',color:['#058DC7', '#50B432', '#e39475', '#DDDF00', '#24CBE5', '#64E572', '#FF9655', '#FFF263', '#6AF9C4'],
label:{normal:{textStyle:{color:'#000000'},formatter: "{b} : {c} ({d}%)"}},radius:'50%',center:['30%','50%'],data:indelexoni2},
{type:'pie',color:['#058DC7', '#50B432', '#e39475', '#DDDF00', '#24CBE5', '#64E572', '#FF9655', '#FFF263', '#6AF9C4'],
label:{normal:{textStyle:{color:'#000000'},formatter: "{b} : {c} ({d}%)"}},radius:'50%',center:['70%','50%'],data:indelfunc2}]
});
indelpieb = document.getElementById("indelpie"),
pindelp=document.getElementById("pindelp"),
nindelp=document.getElementById("nindelp"),
indelpieb.onchange = function() {
		name = indelpieb.options[document.getElementById("indelpie").options.selectedIndex].value, 
		indelfunc2 = indelfuncz2[name], 
		indelexoni2 = indelexoniz2[name], 
		indel_exoni1.setOption({series:[{data:indelexoni2},{data:indelfunc2}]})
	};
	   $("#pindelp").click(function() {
		indelpieb.selectedIndex = indelpieb.selectedIndex - 1,nindelp.style.visibility="visible", -1 == indelpieb.selectedIndex && (indelpieb.selectedIndex = $("#indelpie").find("option:first").index(),pindelp.style.visibility="hidden"), indelpieb.onchange();
	}), $("#nindelp").click(function() {
		indelpieb.selectedIndex = indelpieb.selectedIndex + 1,pindelp.style.visibility="visible", -1 == indelpieb.selectedIndex && (indelpieb.selectedIndex = $("#indelpie").find("option:last").index(),nindelp.style.visibility="hidden"), indelpieb.onchange();
	})};
	if(document.getElementById("INdel_type")){
	indel_type=echarts.init(document.getElementById("INdel_type")),
	indel_type.setOption({
		title: {
			text: "基因组INDEL特征"
		},
		toolbox: {
			show: 1,
			itemSize:30,
			feature: {
				saveAsImage: {
					show: 1
				}
			},
			iconStyle:{emphasis:{borderColor:'#0f5b7c'}}
		},
		tooltip : {
        trigger: 'item',
        formatter: "{b} : {c} ({d}%)"
    },series:[
	{type:'pie',
	label:{normal:{textStyle:{color:'#000000'},formatter: "{b} : {c} ({d}%)"}},
	radius:'50%',
	center : ['30%', '50%'],
	data:indel_homhet},
	{type:'pie',
	label:{normal:{textStyle:{color:'#000000'},formatter: "{b} : {c} ({d}%)"}},
	radius:'50%',
	center : ['70%', '50%'],
	data:indel_newold}
	]
	});
	indelhunb = document.getElementById("indelhun"),
	pinhun=document.getElementById("pinhun"),
	ninhun=document.getElementById("ninhun"),
	indelhunb.onchange = function() {
		name = indelhunb.options[document.getElementById("indelhun").options.selectedIndex].value, 
		indel_homhet = indel_homhetz[name], 
		indel_newold = indel_newoldz[name], 
		indel_type.setOption({
			series:[
	{data:indel_homhet},	{data:indel_newold}]
	})
	};
	
	 $("#pinhun").click(function() {
		indelhunb.selectedIndex = indelhunb.selectedIndex - 1,ninhun.style.visibility="visible", -1 == indelhunb.selectedIndex && (indelhunb.selectedIndex = $("#indelhun").find("option:first").index(),pinhun.style.visibility="hidden"), indelhunb.onchange();
	}), $("#ninhun").click(function() {
		indelhunb.selectedIndex = indelhunb.selectedIndex + 1,pinhun.style.visibility="visible", -1 == indelhunb.selectedIndex && (indelhunb.selectedIndex = $("#indelhun").find("option:last").index(),ninhun.style.visibility="hidden"), indelhunb.onchange();
	})};








	
if(document.getElementById("SV_type")){
	sv_type=echarts.init(document.getElementById("SV_type")),
	sv_type.setOption({
		title: {
			text: "SV检测结果(点击以查看不同区域的SV数量)"
		},
		toolbox: {
			show: 1,
			itemSize:30,
			feature: {
				saveAsImage: {
					show: 1
				}
			},
			iconStyle:{emphasis:{borderColor:'#0f5b7c'}}
		},
		tooltip : {
        trigger: 'item',
        formatter: "{b} : {c} ({d}%)"
    },series:[{type:'pie',label:{normal:{textStyle:{color:'#000000'},formatter: "{b} : {c} ({d}%)"}},radius:'60%',data:sv_data2}]
	}),svpieb = document.getElementById("svpie"),
	psv=document.getElementById("psv"),
	nsv=document.getElementById("nsv"),
	sv_type.on('click',function(param){name = svpieb.options[document.getElementById("svpie").options.selectedIndex].value;if(param.name=='Translocation'){sv_type.setOption({series:[{data:sv_traz[name]}]})}else if(param.name=='Inversion'){sv_type.setOption({series:[{data:sv_invz[name]}]})}else if(param.name=='Deletion'){sv_type.setOption({series:[{data:sv_delz[name]}]})}else if(param.name=='Insertion'){sv_type.setOption({series:[{data:sv_insz[name]}]})}else if(param.name=='ITX'){sv_type.setOption({series:[{data:sv_itx2z[name]}]})}else if(param.name=='CTX'){sv_type.setOption({series:[{data:sv_ctx2z[name]}]})}else if(param.name=='INS'){sv_type.setOption({series:[{data:sv_ins2z[name]}]})}else if(param.name=='INV'){sv_type.setOption({series:[{data:sv_inv2z[name]}]})}else if(param.name=='DEL'){sv_type.setOption({series:[{data:sv_del2z[name]}]})}else if(param.name=='TDP'){sv_type.setOption({series:[{data:sv_tdp2z[name]}]})}else{sv_type.setOption({series:[{data:sv_dataz2[name]}]})}});
	svpieb.onchange = function() {
		name = svpieb.options[document.getElementById("svpie").options.selectedIndex].value, 
		sv_data2 = sv_dataz2[name], 
		sv_type.setOption({series:[{data:sv_data2}]})
	};$("#psv").click(function() {
		svpieb.selectedIndex = svpieb.selectedIndex - 1,nsv.style.visibility="visible", -1 == svpieb.selectedIndex && (svpieb.selectedIndex = $("#svpie").find("option:first").index(),psv.style.visibility="hidden"), svpieb.onchange();
	}), $("#nsv").click(function() {
		svpieb.selectedIndex = svpieb.selectedIndex + 1,psv.style.visibility="visible", -1 == svpieb.selectedIndex && (svpieb.selectedIndex = $("#svpie").find("option:last").index(),nsv.style.visibility="hidden"), svpieb.onchange();
	})};
	if(document.getElementById("CNv_type")){
	cnv_type=echarts.init(document.getElementById("CNv_type")),
	cnv_type.setOption({
		title: {
			text: "CNV检测结果(点击以查看不同区域的CNV数量)"
		},
		toolbox: {
			show: 1,
			itemSize:30,
			feature: {
				saveAsImage: {
					show: 1
				}
			},
			iconStyle:{emphasis:{borderColor:'#0f5b7c'}}
		},
		tooltip : {
        trigger: 'item',
        formatter: "{b} : {c} ({d}%)"
    },series:[{type:'pie',label:{normal:{textStyle:{color:'#000000'},formatter: "{b} : {c} ({d}%)"}},radius:'60%',data:cnv_data2}]
	}),cnvpieb = document.getElementById("cnvpie"),
		pcnv=document.getElementById("pcnv"),
	ncnv=document.getElementById("ncnv"),
	cnv_type.on('click',function(param){name = cnvpieb.options[document.getElementById("cnvpie").options.selectedIndex].value;if(param.name=='gain'){cnv_type.setOption({series:[{data:cnv_gaiz[name]}]})}else if(param.name=='loss'){cnv_type.setOption({series:[{data:cnv_losz[name]}]})}else{cnv_type.setOption({series:[{data:cnv_dataz2[name]}]})}});
	cnvpieb.onchange = function() {
		name = cnvpieb.options[document.getElementById("cnvpie").options.selectedIndex].value, 
		cnv_data2 = cnv_dataz2[name], 
		cnv_type.setOption({series:[{data:cnv_data2}]})
	};$("#pcnv").click(function() {
		cnvpieb.selectedIndex = cnvpieb.selectedIndex - 1,ncnv.style.visibility="visible", -1 == cnvpieb.selectedIndex && (cnvpieb.selectedIndex = $("#cnvpie").find("option:first").index(),pcnv.style.visibility="hidden"), cnvpieb.onchange();
	}), $("#ncnv").click(function() {
		cnvpieb.selectedIndex = cnvpieb.selectedIndex + 1,pcnv.style.visibility="visible", -1 == cnvpieb.selectedIndex && (cnvpieb.selectedIndex = $("#cnvpie").find("option:last").index(),ncnv.style.visibility="hidden"), cnvpieb.onchange();
	})};
	if(document.getElementById("BArplot")){ 
barplot=echarts.init(document.getElementById("BArplot"));
barplot.setOption({
	title: {
			text: '候选基因排序'
		},
		tooltip:{trigger:'item',formatter: "{b} : {c}"},
		toolbox: {
			show: 1,
			itemSize:30,
			feature: {
				saveAsImage: {
					show: 1
				}
			},
			iconStyle:{emphasis:{borderColor:'#0f5b7c'}}
		},yAxis:[{type:'category',inverse:true,data:genebar}],
		visualMap:{min:0,max:1,dimension:'0'}, 
		xAxis:[{type:'value',max:1,min:0,position:'top'}],
		series:[{type:'bar',data:scorebar}]
});
};
setTimeout(function() {
		window.onresize = function() {
			if(document.getElementById("etu")){ute.resize(),f.resize()};
			if(document.getElementById("eta")){ate.resize(),e.resize()};
			if(document.getElementById("etb")){bte.resize(),e.resize()};
			if(document.getElementById("QMh")){qmh.resize()};
			if(document.getElementById("ERm")){erm.resize()};
			if(document.getElementById("GCz")){gcz.resize()};
			if(document.getElementById("ERz")){erz.resize()};
			if(document.getElementById("SNv_exoni1")){snv_exoni1.resize()};
			if(document.getElementById("SNv_type")){snv_type.resize()};
			if(document.getElementById("INdellen")){indellen.resize()};
			if(document.getElementById("INdel_exoni1")){indel_exoni1.resize()};
			if(document.getElementById("INdel_type")){indel_type.resize()};
			if(document.getElementById("SV_type")){sv_type.resize()};
			if(document.getElementById("CNv_type")){cnv_type.resize()};
			if(document.getElementById("BArplot")){barplot.resize()};
		};
	}, 50);
})
