var t = JSON.parse(JSON.parse(document.getElementById("targetid").dataset.targets));



function mir(mi) {
    return `

     <tr>
        <td style="width:220px"><a href="${mi.urlF}" target='_blank'><b>Target Transcript:</b> ${mi.name}</a>&nbsp;&nbsp;&nbsp;&nbsp;</td>
        <td style="width:200px"><b>Range transcript: </b>${mi.Ref_range}&nbsp;&nbsp;&nbsp;&nbsp;</td>
        <td style="width:100px"><b>Score:&nbsp;&nbsp;</b>${mi.score}&nbsp;&nbsp;&nbsp;&nbsp;</td>
        <td style="width:180px"><b>Alignment Length:&nbsp;&nbsp;</b>${mi.aln_len}&nbsp;&nbsp;&nbsp;&nbsp;</td>
        <td style="width:130px"><b>Energy:</b> ${mi.energy}</td>
        <td><a href="${mi.urlm}" target='_blank'><b>MicroRNA:&nbsp;&nbsp;</b> ${mi.mirna}</a>&nbsp;&nbsp;&nbsp;&nbsp;</td>
    </tr>
    <tr><td colspan="6"> <br></td></tr>
    <tr>
    <td colspan="6" align="center">
    <pre><b>${mi.algn}</b></pre>
    </td>
    </tr>
     <tr>
    <td colspan="6">
    <br>
    </td>
    </tr>



      `
}


document.getElementById("mir").innerHTML = `
<table style="float:left">
${t.map(mir).join("")}
</table>

`;


