Shiny.addCustomMessageHandler('updateLog', function(message) {
  // Identifica qual caixa de log será atualizada
  var logBox;
  if (message.target === 'cv_logs') {
    logBox = $('#cv_logs');
  } else if (message.target === 'gei_logs') {
    logBox = $('#gei_logs');
  } else if (message.target === 'predict_logs') {
    logBox = $('#predict_logs');
  } else {
    console.error('Target inválido na mensagem:', message.target);
    return; // Não faz nada se o target for inválido
  }

  if (message.action === 'clearLogs') {
    logBox.text('');
  } else if (message.action === 'addLog') {
    logBox.text(logBox.text() + message.content);
    logBox.scrollTop(logBox[0].scrollHeight);
  } else {
    console.error('Ação inválida na mensagem:', message.action);
  }
});