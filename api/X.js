
// X -- Symplified Scripting Environment

// https://developer.mozilla.org/en-US/docs/Web/JavaScript/Guide/Modules

const visit = window.open;
const cls = console.clear;
const info = console.log;
const warn = console.warn;
const error = console.error;
const group = (o, title) => {
  console.group(title);
  console.log(o);
  console.groupEnd();
}

const X = 
{
logo: '𝕏',
utf: 'U+1D52B',
about: '𝕏 -- Scripting API, Version 2024-MAR-16',
home: 'http://dave-ryzen/omega/api',
dunsel: ()=>{}
, //---------------------------
todo: (method)=>{
    warn(`TODO: 𝕏.${method}`);
    return X;
}
, //---------------------------
help: ()=>{
    visit(`${X.home}/help/X`);
    return X;
}
, //---------------------------
toc: ()=>{
    visit(`${X.home}/help/X/toc.html`);
    return X;
}
, //---------------------------
init: ()=>{ // 𝕏 U+1D52B 
  info("𝕏 Initialized");
  return X;
}
, //---------------------------
list: ()=>{
    const scripts = [ ... document.scripts ];
    const descriptors = scripts.map(s=>{
        return { 
            url: s.src,
            embedded: X.isEmbedded(s),
            system: X.isSystem(s)
        };        
    });
    group(descriptors, "Scripts");
    return descriptors;
}
, //---------------------------
count: ()=>{
    return document.scripts.length;
}
, //---------------------------
select: (regex)=>{
    const scripts = [ ... document.scripts ];
    return scripts.filter(o=>o.src.match(regex));
}
, //---------------------------
find: (url)=>{
    const scripts = [ ... document.scripts ];
    const match = scripts.filter(o=>o.src==url);
    return match[0];
}
, //---------------------------
load: (url, callback)=>{
    callback = callback || (()=>{});
    let script = X.find(url);
    if (script) {
        callback(script);
    } else {
        script = document.createElement("script");
        script.onload = callback;
        document.body.appendChild(script);
        script.src = url;
    }
    return X;
}
, //---------------------------
unload: (url)=>{
    let script = X.find(url);
    if (script) {
        if (X.isSystem(script)) {
            warn(`System script can't be unloaded:`, { url });
        } else {
            const parent = script.parentElement;
            if (parent) parent.removeChild(script);
            info(`Script unloaded:`, { url });    
        }
    } else {
        warn(`Script not found:`, { url });
    }
    return X;
}
, //---------------------------
isSystem: (script) => {
    return script.hasAttribute('system');
}
, //---------------------------
isEmbedded: (script) => {
    return script.src.length ? false : true;
}
, //---------------------------
isType: (script, type) => {
    return script.getAttribute('type') == type;
}
, //---------------------------
isModule: (script) => {
    return this.isType(script, "module");
}
, //---------------------------
isImportMap: (script) => {
    return this.isType(script, "importmap");
}
};

window.addEventListener('load', (evt)=>{
    X.init();
});

