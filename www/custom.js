Shiny.addCustomMessageHandler('updateLog', function(message) {
  var logBox = $('#cv_logs');
  if (message === 'clearLogs') {
    logBox.text('');
  } else {
    logBox.text(logBox.text() + message);
    logBox.scrollTop(logBox[0].scrollHeight); // Auto-scroll to the bottom
  }
});

Shiny.addCustomMessageHandler('updateLog', function(message) {
  var logBox2 = $('#gei_logs');
  if (message === 'clearLogs') {
    logBox2.text('');
  } else {
    logBox2.text(logBox2.text() + message);
    logBox2.scrollTop(logBox2[0].scrollHeight); // Auto-scroll to the bottom
  }
});