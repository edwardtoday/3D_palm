mw.loader.implement("ext.gadget.site-lib",function($){var wgUXS=function(wg,hans,hant,cn,tw,hk,sg,zh,mo,my){var ret={'zh':zh||hans||hant||cn||tw||hk||sg||mo||my,'zh-hans':hans||cn||sg||my,'zh-hant':hant||tw||hk||mo,'zh-cn':cn||hans||sg||my,'zh-sg':sg||hans||cn||my,'zh-tw':tw||hant||hk||mo,'zh-hk':hk||hant||tw||mo}
return ret[wg]||zh||hans||hant||cn||tw||hk||sg||mo||my;}
window.wgULS=function(hans,hant,cn,tw,hk,sg,zh,mo,my){return wgUXS(wgUserLanguage,hans,hant,cn,tw,hk,sg,zh,mo,my);};window.wgUVS=function(hans,hant,cn,tw,hk,sg,zh,mo,my){return wgUXS(wgUserVariant,hans,hant,cn,tw,hk,sg,zh,mo,my);};window.getParamValue=mw.util.getParamValue;;},{},{});

/* cache key: zhwiki:resourceloader:filter:minify-js:4:6a3ee7a010c75302fd500091e82f1889 */