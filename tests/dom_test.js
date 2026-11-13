const { JSDOM } = require('jsdom');
const { spawn } = require('child_process');
const P = 8101, BASE = 'http://127.0.0.1:' + P + '/';
let fails = 0;
const chk = (c, m) => { console.log('  ' + (c ? 'pass' : 'FAIL') + '  ' + m); if (!c) fails++; };
const sleep = ms => new Promise(r => setTimeout(r, ms));

(async () => {
  const srv = spawn('./crispr', ['serve', String(P)], { cwd: process.env.CRISPR_DIR || process.cwd() });
  await sleep(600);
  try {
    const dom = await JSDOM.fromURL(BASE, {
      runScripts: 'dangerously', pretendToBeVisual: true,
      beforeParse(w) {
        w.Element.prototype.scrollIntoView = () => {};
        w.fetch = async (u, o = {}) => {
          let body = o.body;
          if (body && typeof body !== 'string') body = await new Promise(res => { const fr = new w.FileReader(); fr.onload = () => res(Buffer.from(fr.result)); fr.readAsArrayBuffer(body); });
          return fetch(new URL(u, BASE), { method: o.method, body });
        };
      }
    });
    const w = dom.window, d = w.document, $ = id => d.getElementById(id);
    const until = async (f, m, t = 90000) => { const t0 = Date.now(); while (Date.now() - t0 < t) { try { if (f()) return true; } catch (e) {} await sleep(50); } chk(false, 'timeout: ' + m); return false; };
    const click = id => $(id).dispatchEvent(new w.MouseEvent('click', { bubbles: true }));
    const err = n => $('err' + n).textContent;

    await until(() => $('status').textContent.includes('None loaded'), 'initial status');
    chk(d.querySelectorAll('#nav button').length === 8 && d.querySelectorAll('section').length === 8, '8 sidebar pages and 8 sections');
    chk(d.querySelectorAll('details').length === 8, 'every page has a collapsible "How it works" note');
    d.querySelectorAll('nav button').forEach((b, i) => { b.click(); chk($('p' + (i + 1)).classList.contains('on'), 'nav to page ' + (i + 1)); });

    // error paths with no genome
    for (const [btn, n] of [['b2', 2], ['b3', 3], ['b4', 4], ['b5', 5], ['b7', 7]]) {
      click(btn); await until(() => err(n).length > 0, 'error page ' + n);
      chk(/No genome is loaded/.test(err(n)), 'page ' + n + ' with no genome shows plain error text');
    }
    chk([...d.styleSheets[0].cssRules].some(r => r.cssText.startsWith('.err') && r.cssText.includes('var(--warn)')) && [...d.styleSheets[0].cssRules].find(r => r.cssText.startsWith(':root')).cssText.includes('--warn: #a8372c'), 'error text uses the muted red (checked in the stylesheet; jsdom cannot resolve var())');

    // page 1
    const synth = d.querySelector('#p1 button[data-n="300000"]');
    synth.click();
    chk(synth.disabled === true, 'button disabled while the request runs');
    await until(() => $('out1').textContent.includes('Loaded'), 'synthetic load');
    chk(synth.disabled === false, 'button enabled again afterwards');
    chk($('status').textContent.includes('300,000') && $('from').value === '85000' && $('to').value === '115000', 'sidebar status filled, region prefilled around planted guide');

    // page 2
    click('b2'); await until(() => $('out2').textContent.includes('All three agree'), 'page 2');
    chk($('out2').querySelectorAll('table')[0].querySelectorAll('tr').length === 4 && $('out2').querySelectorAll('svg').length === 2, 'page 2: 3 matcher rows, two bar charts (2-letter and 16-letter)');

    // page 3
    click('b3'); await until(() => $('out3').querySelectorAll('tr.pick').length === 13, 'page 3');
    chk($('out3').querySelectorAll('tr.pick').length === 13, 'page 3: 10 best + 3 worst rows, strand map drawn: ' + ($('out3').querySelectorAll('svg').length === 1));
    const planted = [...$('out3').querySelectorAll('tr.pick')].slice(10).find(tr => tr.textContent.includes('100,000'));
    chk(!!planted, 'planted guide (position 100,000) is in the worst table');
    planted.dispatchEvent(new w.MouseEvent('click', { bubbles: true }));
    await until(() => $('hits3').textContent.includes('13 sites'), 'hits panel');
    chk($('hits3').querySelectorAll('tbody tr, table tr').length - 1 === 13 && $('hits3').querySelectorAll('svg rect').length >= 20, 'hits: 13 rows, mismatch chart with 20 bars');
    chk($('hits3').querySelectorAll('.mm').length > 0, 'mismatched letters are highlighted');
    $('to').value = '10'; $('from').value = '500'; click('b3');
    await until(() => /23 bases/.test(err(3)), 'region error');
    chk(/23 bases/.test(err(3)), 'page 3: bad region shows plain error: ' + err(3).slice(0, 55));
    $('from').value = '85000'; $('to').value = '115000';

    // page 4
    click('b4'); await until(() => $('out4').textContent.includes('Identical results: 10 of 10'), 'page 4');
    chk(/faster per query/.test($('out4').textContent) && $('out4').querySelectorAll('tr').length === 11, 'page 4: speed-up text, 10 guide rows');
    console.log('       ' + $('out4').querySelector('p').textContent.slice(0, 230));

    // page 5
    click('b5'); await until(() => $('out5').textContent.includes('islands found'), 'page 5');
    chk($('out5').querySelectorAll('tr').length === 4 && $('out5').querySelectorAll('svg').length === 1, 'page 5: GC track plus 3 island rows');

    // page 6
    click('b6'); await until(() => $('out6').querySelectorAll('.ablock').length > 0, 'page 6');
    chk(/Same optimal score/.test($('out6').textContent) && $('out6').querySelectorAll('svg').length === 1, 'page 6: scores equal, memory bars, alignment blocks');
    $('a6').value = 'GATTAXA'; $('b6t').value = 'GCATGCA'; click('b6c'); await until(() => err(6).length > 0, 'page 6 bad input');
    chk(/A, C, G or T/.test(err(6)), 'page 6: bad letters give plain error');
    $('a6').value = 'GATTACA'; click('b6c'); await until(() => /Same optimal score/.test($('out6').textContent) && err(6) === '', 'page 6 custom');
    // (empty b sequence is rejected too)
    chk(err(6) === '', 'page 6: error text cleared on the next successful run');

    // page 7
    click('b7'); await until(() => $('out7').querySelectorAll('.tag.pass').length === 3, 'page 7 vault');
    chk($('out7').querySelectorAll('.flow > div').length === 4 && $('out7').querySelectorAll('.tag.fail').length === 0, 'page 7: four-step flow, 3 PASS, 0 FAIL');
    click('b7a'); await until(() => $('att7').textContent.includes('What this does not prove'), 'page 7 attack');
    chk($('att7').querySelectorAll('svg').length === 1 && $('att7').querySelectorAll('tr').length === 12, 'page 7: rho chart and 11 table rows');
    chk($('att7').querySelectorAll('polyline').length === 2, 'page 7: measured curve and theory line drawn');

    // page 8
    click('b8'); await until(() => $('out8').querySelectorAll('svg').length === 4, 'page 8');
    chk($('out8').querySelectorAll('svg').length === 4 && $('out8').querySelectorAll('table').length === 2, 'page 8: 4 charts, 2 tables');

    // garbage upload through the page's own upload code path
    await w.upload(new w.Blob(['this is not a genome'], { type: 'text/plain' }));
    chk(/FASTA/.test(err(1)), 'page 1: garbage upload shows plain error: ' + err(1).slice(0, 60));
    const okfa = '>t\nACGT'.padEnd(2000, 'ACGTTGCA');
    await w.upload(new w.Blob([okfa]));
    chk(err(1) === '' && /Loaded/.test($('out1').textContent), 'page 1: valid upload loads (small genome)');

    // server stopped: page should say so plainly
    srv.kill(); await sleep(400);
    click('b2'); await until(() => /Cannot reach the server/.test(err(2)), 'server down');
    chk(/Cannot reach the server/.test(err(2)), 'server down: plain message');
  } catch (e) { console.log('TEST CRASH', e); fails++; }
  srv.kill();
  console.log('\nDOM test:', fails ? fails + ' FAILED' : 'ALL PASSED');
  process.exit(fails ? 1 : 0);
})();
